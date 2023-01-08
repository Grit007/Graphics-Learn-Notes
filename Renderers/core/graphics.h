#ifndef RENDERERS_GRAPHI_H
#define RENDERERS_GRAPHI_H

#include "../shader/shader.h"

Matrix viewport(int x, int y, int w, int h) {
    Viewport = Matrix::identity();
    Viewport[0][3] = x+w/2.f;
    Viewport[1][3] = y+h/2.f;
    Viewport[2][3] = 1.f;
    Viewport[0][0] = w/2.f;
    Viewport[1][1] = h/2.f;
    Viewport[2][2] = 0;
    return Viewport;
}

Matrix projection(float coeff) {
    Projection = Matrix::identity();
    Projection[3][2] = coeff;
    return Projection;
}

Matrix lookat(Vec3f eye, Vec3f target, Vec3f up) {
    Vec3f z = (eye-target).normalize();
    Vec3f x = cross(up,z).normalize();
    Vec3f y = cross(z,x).normalize();
    Matrix Minv = Matrix::identity();
    Matrix Tr   = Matrix::identity();
    for (int i=0; i<3; i++) {
        Minv[0][i] = x[i];
        Minv[1][i] = y[i];
        Minv[2][i] = z[i];
        Tr[i][3] = -target[i];
    }
    ModelView = Minv*Tr;
    return ModelView;
}

Vec3f barycentric(Vec2f A, Vec2f B, Vec2f C, Vec2f P) {
    Vec3f s[2];
    for (int i=2; i--; ) {
        s[i][0] = C[i]-A[i];
        s[i][1] = B[i]-A[i];
        s[i][2] = A[i]-P[i];
    }
    Vec3f u = cross(s[0], s[1]);
    if (std::abs(u[2])>1e-2) // dont forget that u[2] is integer. If it is zero then triangle ABC is degenerate
        return Vec3f(1.f-(u.x+u.y)/u.z, u.y/u.z, u.x/u.z);
    return Vec3f(-1,1,1); // in this case generate negative coordinates, it will be thrown away by the rasterizator
}

void DrawPoint(int x1, int y1, TGAImage &image, TGAColor color) {
    image.set(x1, y1, color);
}

void DrawLine(int x1, int y1, int x2, int y2, TGAImage &image, TGAColor color){
    float dt = 0.01;
    int dx = x2 - x1;
    int dy = y2 - y1;
    for(float t = 0.0; t <= 1.0; t += dt){
        int x = x1 + dx * t;
        int y = y1 + dy * t;
        image.set(x, y, color);
    }
}

void DrawLineDDA(int x1, int y1, int x2, int y2, TGAImage &image, TGAColor color){
    int dx = x2 - x1;
    int dy = y2 - y1;
    int steps, k;
    float xIncrement, yIncrement, x = x1, y = y1;
    if(fabs(dx) > fabs(dy))
        steps = fabs(dx);
    else
        steps = fabs(dy);
    xIncrement = float(dx) / float(steps);
    yIncrement = float(dy) / float(steps);
    image.set(x, y, color);
    for(k = 0; k < steps; k++){
        x += xIncrement;
        y += yIncrement;
        image.set(x, y, color);
    }
}

void DrawLineBresenham(int x1, int y1, int x2, int y2, TGAImage &image, TGAColor color){
    bool steep = false;
    // 判断直线是否"陡峭"
    if(abs(x1-x2) < abs(y1-y2)){
        std::swap(x1, y1);
        std::swap(x2,y2);
        steep = true;
    }
    // 确保直线起始点小于结束点
    if(x1 > x2){
        std::swap(x1, x2);
        std::swap(y1,y2);
    }
    int y = y1, eps = 0, yi = 1;
    int dx = x2 - x1;
    int dy = y2 - y1;
    // 处理 [-1, 0] 范围内的斜率
    if (dy < 0) {
        yi = -1;
        dy = -dy;
    }
    if(steep){
        for(int x = x1; x <= x2; x++){
            image.set(y, x, color);
            eps += dy;
            if((eps << 1) >= dx)  {
                y = y + yi;
                eps -= dx;
            }
        }
    } else {
        for(int x = x1; x <= x2; x++){
            image.set(x, y, color);
            eps += dy;
            if((eps << 1) >= dx)  {
                y = y + yi;
                eps -= dx;
            }
        }
    }
}

void DrawWireframe(Vec3f *pts, TGAImage &image, TGAColor color){
    for (int j=0; j<3; j++){
        Vec3f p1 = pts[j];
        Vec3f p2 = pts[(j+1)%3];
        DrawLineBresenham(p1.x, p1.y, p2.x, p2.y, image, color);
    }
}

void DrawTriangle(mat<4,3,float> &clipc, IShader &shader, TGAImage &image, float *zbuffer, Camera &camera) {
    mat<3,4,float> pts  = (Viewport*clipc).transpose(); // transposed to ease access to each of the points
    mat<3,2,float> pts2;
    for (int i=0; i<3; i++) pts2[i] = proj<2>(pts[i]/pts[i][3]);
    // 求三角形包围盒
    Vec2f bboxmin( std::numeric_limits<float>::max(),  std::numeric_limits<float>::max());
    Vec2f bboxmax(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());
    Vec2f clamp(image.get_width()-1, image.get_height()-1);
    for (int i=0; i<3; i++) {
        for (int j=0; j<2; j++) {
            bboxmin[j] = std::max(0.f,      std::min(bboxmin[j], pts2[i][j]));
            bboxmax[j] = std::min(clamp[j], std::max(bboxmax[j], pts2[i][j]));
        }
    }
    Vec2i P;
    TGAColor color;
    // 光栅化
    for (P.x=bboxmin.x; P.x<=bboxmax.x; P.x++) {
        for (P.y=bboxmin.y; P.y<=bboxmax.y; P.y++) {
            Vec3f bc_screen  = barycentric(pts2[0], pts2[1], pts2[2], P);
            // 透视校正
            Vec3f bc_clip    = Vec3f(bc_screen.x/pts[0][3], bc_screen.y/pts[1][3], bc_screen.z/pts[2][3]);
            bc_clip = bc_clip/(bc_clip.x+bc_clip.y+bc_clip.z);
            float frag_depth = clipc[2]*bc_clip;
            // 深度测试
            if (bc_screen.x<0 || bc_screen.y<0 || bc_screen.z<0 || zbuffer[P.x+P.y*image.get_width()]>frag_depth) continue;
            // 着色
            bool discard = shader.fragment(bc_clip, color, camera);
            if (!discard) {
                zbuffer[P.x+P.y*image.get_width()] = frag_depth;
                image.set(P.x, P.y, color);
            }
        }
    }
}

#endif //RENDERERS_GRAPHI_H
