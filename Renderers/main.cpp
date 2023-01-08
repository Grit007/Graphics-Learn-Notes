#include <limits>
#include "platform/win32.h"
#include "core/camera.h"
#include "core/graphics.h"

void ClearZbuffer(int width, int height, float* zbuffer, float* shadowbuffer)
{
    for (int i=width*height; --i; ) {
        zbuffer[i] = shadowbuffer[i] = -std::numeric_limits<float>::max();
    }
}

void ClearFramebuffer(int width, int height, unsigned char* framebuffer)
{
    memset(framebuffer, 50, width * height * 4);
}

void UpdateMatrix(Camera &camera){
    ModelView = lookat(camera.eye, camera.target, camera.up);
    Projection = projection(-1.f/(camera.eye-camera.target).norm());
}

int main(int argc, char** argv) {
//    model = new Model("../obj/african_head/african_head.obj");
    model = new Model("../obj/pbr/backpack.obj");
    // create camera
    Camera camera(eye, target, up, (float)(width) / height);
    // init MVP matrix
    ModelView = lookat(eye, target, up);
    Projection = projection(-1.f/(eye-target).norm());
    Viewport = viewport(width/8, height/8, width*3/4, height*3/4);
    // transform light_dir to clip space
    light_dir = proj<3>(Projection*ModelView*embed<4>(light_dir, 0.f));
    // initialize window
    window_init(width, height, "Renderers");
    // set draw mode
    DrawMode dm = Triangle; // {Point, Line, WireFrame, Triangle}
    // init shader type
//    BlinnPhongShader shader;
    DepthShader depthShader;
    PBRShader shader;
    light_dir_depth.normalize();
    Matrix depthM, normalM;
    std::vector<int> face;
    Vec3f pts[3];
    // render loop
    while(!window->is_close){
        depth.clear();
        // clear buffer
        ClearZbuffer(width, height, zbuffer, shadowbuffer);
        ClearFramebuffer(width, height, image.buffer());
        // handle events and update view, perspective matrix
        handle_events(camera);
        UpdateMatrix(camera);
        switch (dm) {
            case Point:
                DrawPoint(200, 200, image, red);
                break;
            case Line:
//                DrawLine(80, 40, 200, 200, image, red);
//                DrawLineDDA(80, 40, 200, 200, image, red);
                DrawLineBresenham(80, 40, 200, 200, image, red);
                break;
            case WireFrame:
            {
//                    for(int m=1; m<argc; m++) {
//                        model = new Model(argv[m]);
                for (int i = 0; i < model->nfaces(); i++) {
                    face = model->face(i);
                    for (int j = 0; j < 3; j++) {
                        pts[j] = proj<3>(Viewport * Projection * ModelView *
                                         embed<4>(model->vert(face[j])));
                    }
                    DrawWireframe(pts, image, red);
                }
//                    }
            }
                break;
            case Triangle:
            {
//                    for(int m=1; m<argc; m++){
//                        model = new Model(argv[m]);
                {
                    // rendering the shadow buffer
                    lookat(light_dir_depth, camera.target, camera.up);
                    projection(0);
                    viewport(width/8, height/8, width*3/4, height*3/4);

                    for (int i=0; i<model->nfaces(); i++) {
                        for (int j=0; j<3; j++) {
                            depthShader.vertex(i, j);
                        }
                        DrawTriangle(depthShader.varying_tri, depthShader, depth, shadowbuffer, camera);
                    }
                }
                depthM = Projection*ModelView;
                { // rendering the frame buffer
                    lookat(camera.eye, camera.target, camera.up);
                    projection(-1.f/((camera.eye-camera.target).norm()));
                    viewport(width/8, height/8, width*3/4, height*3/4);

                    Mshadow = depthM*(Projection*ModelView).invert();
                    normalM = (Projection*ModelView).invert_transpose();
//                    transform matrix 4x4 to 3x3 ;
                    for(int i=0; i<3; i++)
                        for(int j=0; j<3; j++)
                            NormalMatrix[i][j] = normalM[i][j];

                    for (int i=0; i<model->nfaces(); i++) {
                        for (int j=0; j<3; j++) {
                            shader.vertex(i, j);
                        }
                        DrawTriangle(shader.varying_tri, shader, image, zbuffer, camera);
                    }
                }
//                    }
            }
                break;
        }
        // reset mouse information
        window->mouse_info.wheel_delta = 0;
        window->mouse_info.orbit_delta = Vec2f(0,0);
        window->mouse_info.fv_delta = Vec2f(0, 0);
//        to have the origin at the left bottom corner of the image
        depth.flip_vertically();
        depth.write_tga_file("../depth.tga");
        image.flip_vertically();
        image.write_tga_file("../output.tga");
        // send framebuffer to window to display
        window_draw(image.buffer());
        msg_dispatch();
    }

    image.write_tga_file("../output.tga");
    window_destroy();
    delete model;
    delete [] zbuffer;
    delete [] shadowbuffer;
    return 0;
}