#ifndef RENDERERS_SHADER_H
#define RENDERERS_SHADER_H

#include "../common.h"

class IShader {
public:
    mat<2,3,float> varying_uv;  // triangle uv coordinates, written by the vertex shader, read by the fragment shader
    mat<4,3,float> varying_tri; // triangle coordinates (clip coordinates), written by VS, read by FS
    mat<3,3,float> varying_nrm; // normal per vertex to be interpolated by FS
    mat<3,3,float> ndc_tri;     // triangle in normalized device coordinates
    mat<3,3,float> fs_pos; // clip coordinates 3X3 计算半程向量使用

    virtual ~IShader();
    virtual Vec4f vertex(int iface, int nthvert) { return Vec4f();};
    virtual bool fragment(Vec3f bar, TGAColor &color, Camera &camera) { return false;};
};

IShader::~IShader() {}

mat<3,3,float> TBN(mat<3,3,float> ndc_tri1, mat<2,3,float> varying_uv1, Vec3f bn1){
    // 计算TBN矩阵
    mat<3,3,float> A;
    A[0] = ndc_tri1.col(1) - ndc_tri1.col(0);
    A[1] = ndc_tri1.col(2) - ndc_tri1.col(0);
    A[2] = bn1;
    mat<3,3,float> AI = A.invert();
    Vec3f i = AI * Vec3f(varying_uv1[0][1] - varying_uv1[0][0], varying_uv1[0][2] - varying_uv1[0][0], 0);
    Vec3f j = AI * Vec3f(varying_uv1[1][1] - varying_uv1[1][0], varying_uv1[1][2] - varying_uv1[1][0], 0);
    mat<3,3,float> B;
    B.set_col(0, i.normalize());
    B.set_col(1, j.normalize());
    B.set_col(2, bn1);
    return B;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
class DepthShader : public IShader {
public:
    virtual Vec4f vertex(int iface, int nthvert) {
        Vec4f gl_Vertex = Projection*ModelView*embed<4>(model->vert(iface, nthvert));
        varying_tri.set_col(nthvert, gl_Vertex);
        ndc_tri.set_col(nthvert, proj<3>(gl_Vertex/gl_Vertex[3]));
        return gl_Vertex;
    }
    virtual bool fragment(Vec3f bar, TGAColor &color, Camera &camera) {
        color = TGAColor(255, 255, 255)*((ndc_tri*bar).z);
        return false;
    }
};

class BlinnPhongShader : public IShader {
public:
    virtual Vec4f vertex(int iface, int nthvert) {
        varying_uv.set_col(nthvert, model->uv(iface, nthvert));
        varying_nrm.set_col(nthvert, proj<3>((Projection*ModelView).invert_transpose()*embed<4>(model->normal(iface, nthvert), 0.f)));
        Vec4f gl_Vertex = Projection*ModelView*embed<4>(model->vert(iface, nthvert));
        varying_tri.set_col(nthvert, gl_Vertex);
        // 开启布林冯着色时保存 fs_pos
        if(blinn) fs_pos.set_col(nthvert, proj<3>(gl_Vertex));
        ndc_tri.set_col(nthvert, proj<3>(gl_Vertex/gl_Vertex[3]));
        return gl_Vertex;
    }
    virtual bool fragment(Vec3f bar, TGAColor &color, Camera &camera) {
        Vec4f sb_p = Mshadow*embed<4>(varying_tri*bar); // corresponding point in the shadow buffer
        sb_p = sb_p/sb_p[3];
        sb_p =  Viewport*sb_p;
        int idx = int(sb_p[0]) + int(sb_p[1])*width; // index in the shadowbuffer array
        float shadow = .3+.7*(shadowbuffer[idx]<sb_p[2]); // magic coeff to avoid z-fighting
        Vec3f bn = (varying_nrm*bar).normalize();
        Vec2f uv = varying_uv*bar;
        Vec3f curPos = fs_pos * bar;
        // ndc坐标的z值即为深度 现在的范围是-1-1，不缩放到0-1也没事，因为并不是从shadowmap中采样的深度 如果用shadowmap采样则需要缩放
        Vec3f viewDir = (proj<3>(Projection*ModelView*embed<4>(camera.eye)) - curPos).normalize();
        Vec3f lightDir = (light_dir - curPos).normalize();
        mat<3,3,float> TBNMatrix = TBN(ndc_tri, varying_uv, bn);
        Vec3f n = (NormalMatrix*TBNMatrix*model->normal(uv)).normalize();
        // diffuse intensity
        float diff = std::max(0.f, n*lightDir);
        // diffuse color
        color = model->diffuse(uv);
        for(int i=0; i<3; i++) color[i] = pow(color[i]/255.0, gamma)*255.0;
        // reflected light direction
        Vec3f r = (n*(n*lightDir*2.f)-lightDir).normalize();
        float spec = 0.f;
        if(blinn){
            Vec3f halfwayDir = (lightDir + viewDir).normalize();
            /* 半程向量与表面法线的夹角通常会小于观察与反射向量的夹角。所以，如果你想获得和冯氏着色类似的效果，
             * 就必须在使用Blinn-Phong模型时将镜面反光度设置更高一点。通常我们会选择冯氏着色时反光度分量的2到4倍。*/
            spec = pow(std::fmax(n*halfwayDir, 0.0), 4*model->specular(uv));
        }
        else{
            spec = pow(std::fmax(viewDir*r, 0.0), model->specular(uv));
        }
//       ambient + diff + spec, clamp the result and gamma correcting
        for (int i=0; i<3; i++)
        {
            color[i] = pow(std::min<float>(3. + color[i]*shadow*(1.2*diff + .8*spec), 255)/255.0, 1.0/gamma) * 255;
        }
        return false;
    }
};
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// 几何函数(G项) 描述的是微平面自遮挡的情况(自阴影) 若表面被遮挡 则会减少该表面所反射的光线
float GeometrySchlickGGX(float NdotV, float roughness)
{
    float r = (roughness + 1.0);
    float k = (r*r) / 8.0;

    float nom   = NdotV;
    float denom = NdotV * (1.0 - k) + k;

    return nom / denom;
}
float GeometrySmith(Vec3f N, Vec3f V, Vec3f L, float roughness)
{
    float NdotV = std::fmax(N*V, 0.0);
    float NdotL = std::fmax(N*L, 0.0);
    float ggx2 = GeometrySchlickGGX(NdotV, roughness);
    float ggx1 = GeometrySchlickGGX(NdotL, roughness);

    return ggx1 * ggx2;
}
// 法线分布函数(D项) 根据粗糙度来估算微平面的情况 当H趋向一致时 NDF会形成一个明亮的斑点; 当H较为分散时 表现效果看起来就会偏灰暗
float DistributionGGX(Vec3f N, Vec3f H, float roughness)
{
    float a = roughness*roughness;
    float a2 = a*a;
    float NdotH = std::fmax(N*H, 0.0);
    float NdotH2 = NdotH*NdotH;

    float nom   = a2;
    float denom = (NdotH2 * (a2 - 1.0) + 1.0);
    denom = PI * denom * denom;

    return nom / denom;
}
// 菲涅尔方程(F项)计算出的系数其实就是kS 即镜面反射占比部分
Vec3f fresnelSchlick(float cosTheta, Vec3f F0)
{
    return F0 - (F0 - 1.0f) * pow(Clamp(1.0f - cosTheta, 0.0f, 1.0f), 5.0f);
    // 因为运算符重载的缘故，换成上面的写法
//    return F0 + (1.0 - F0) * pow(Clamp(1.0 - cosTheta, 0.0, 1.0), 5.0);
}

Vec3f color2Vec(TGAColor color){
    return Vec3f(pow(color.bgra[0]/255.0, gamma)*255,pow(color.bgra[1]/255.0, gamma)*255,pow(color.bgra[2]/255.0, gamma)*255);
}

float ACESToneMapping(float color, float adapted_lum)
{
    const float A = 2.51f;
    const float B = 0.03f;
    const float C = 2.43f;
    const float D = 0.59f;
    const float E = 0.14f;

    color *= adapted_lum;
    return (color * (A * color + B)) / (color * (C * color + D) + E);
}

class PBRShader : public IShader {
public:
    virtual Vec4f vertex(int iface, int nthvert) {
        varying_uv.set_col(nthvert, model->uv(iface, nthvert));
        varying_nrm.set_col(nthvert, proj<3>((Projection*ModelView).invert_transpose()*embed<4>(model->normal(iface, nthvert), 0.f)));
        Vec4f gl_Vertex = Projection*ModelView*embed<4>(model->vert(iface, nthvert));
        varying_tri.set_col(nthvert, gl_Vertex);
        ndc_tri.set_col(nthvert, proj<3>(gl_Vertex/gl_Vertex[3]));
        return gl_Vertex;
    }
    virtual bool fragment(Vec3f bar, TGAColor &color, Camera &camera) {
        Vec3f bn = (varying_nrm*bar).normalize();
        Vec2f uv = varying_uv*bar;
        Vec3f curPos = fs_pos * bar;
        Vec3f viewDir = (proj<3>(Projection*ModelView*embed<4>(camera.eye)) - curPos).normalize();
        mat<3,3,float> TBNMatrix = TBN(ndc_tri, varying_uv, bn);
        Vec3f n = (NormalMatrix*TBNMatrix*model->normal(uv)).normalize();
        // calculate reflectance at normal incidence; if dia-electric (like plastic) use F0 of 0.04 and if it's a metal, use the albedo color as F0 (metallic workflow)
        Vec3f F0 = Vec3f(0.04f, 0.04f, 0.04f);
        Vec3f albedo = color2Vec(model->diffuse(uv)); // color2Vec这个函数已经做了gamma2line空间的转换了
        float metallic = model->metallic(uv);
        F0 = Mix(F0, albedo, metallic);
        // reflectance equation
        Vec3f Lo = Vec3f(0.0, 0.0, 0.0);
        for(int i=0; i<4; i++)
        {
            Vec3f lightDir = (light_dirs[i] - curPos).normalize();
            Vec3f halfwayDir = (lightDir + viewDir).normalize();
            float distance = Length(light_dirs[i] - curPos);
            float attenuation = 1.0 / (distance * distance);
            Vec3f radiance = light_colors[i] * attenuation;
            // Cook-Torrance BRDF
            float roughness = model->roughness(uv);
            float NDF = DistributionGGX(n, halfwayDir, roughness);
            float G   = GeometrySmith(n, viewDir, lightDir, roughness);
            Vec3f F    = fresnelSchlick(Clamp(halfwayDir*viewDir, 0.0, 1.0), F0);
            Vec3f numerator    = F * NDF * G;
            float denominator = 4.0 * std::fmax(n*viewDir, 0.0) * std::fmax(n*lightDir, 0.0) + 0.0001; // + 0.0001 to prevent divide by zero
            Vec3f specular = numerator / denominator;
            // kS is equal to Fresnel
            Vec3f kS = F;
            Vec3f kD = Vec3f(1.0, 1.0, 1.0) - kS;
            kD = kD * (1.0 - metallic);
            // scale light by NdotL
            float NdotL = std::fmax(n*lightDir, 0.0);
            Lo += Vecmul((Vecmul(kD, albedo) / PI + specular), radiance) * NdotL;
        }
        float ao = model->ao(uv);
        Vec3f ambient = Vec3f(albedo[0]*0.03, albedo[1]*0.03, albedo[2]*0.03) * ao;
        for(int i=0; i<3; i++)
        {
            color[i] = ambient[i] + Lo[i];
            // gamma correct and tonemapping
            color[i] = pow(ACESToneMapping(color[i]/255.0, 0.6), 1.0/gamma)*255;

        }
        return false;
    }
};

#endif //RENDERERS_SHADER_H
