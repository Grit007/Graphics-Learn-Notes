#ifndef RENDERERS_COMMON_H
#define RENDERERS_COMMON_H
#include "utils/maths.h"
#include "utils/model.h"

// 存放公共变量的定义
Matrix ModelView;
Matrix Projection;
Matrix Viewport;
mat<3,3,float> NormalMatrix; // 法线变换矩阵
Matrix Mshadow; // transform  coordinates to shadowbuffer space

Model *model = NULL;
TGAColor red = TGAColor(255,0,0,255);

const int width  = 400;
const int height = 400;

TGAImage image(width, height, TGAImage::RGBA);
TGAImage depth(width, height, TGAImage::RGBA);

const float gamma = 2.2;
float *shadowbuffer   = new float[width*height];
float *zbuffer = new float[width * height];
bool blinn = false; // 是否开启布林冯着色 默认不开启 即采用冯氏着色 按B键可切换
// 布林冯
Vec3f light_dir(1,1,1); // light source(position)
Vec3f light_dir_depth(1,1,1); // light source(position) 阴影pass的时候使用
// pbr
Vec3f light_dirs[4] = {Vec3f(1,1,1), Vec3f(1,-1,-1), Vec3f(-1,1,1), Vec3f(-1,-1,-1)};
Vec3f light_colors[4] = {Vec3f(10,10,10), Vec3f(10,10,10), Vec3f(10,10,10), Vec3f(10,10,10)};

Vec3f       eye(1,1,3); // camera position
Vec3f    target(0,0,0); // camera direction
Vec3f        up(0,1,0); // camera up vector

enum DrawMode {Point, Line, WireFrame, Triangle};

#endif //RENDERERS_COMMON_H
