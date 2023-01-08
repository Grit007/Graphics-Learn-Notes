#ifndef RENDERERS_CAMERA_H
#define RENDERERS_CAMERA_H

#include "../utils/maths.h"

#define PI 3.1415926
#define EPSILON 1e-5f

class Camera
{
public:
    Camera(Vec3f e, Vec3f t, Vec3f up, float aspect);
    ~Camera();

    Vec3f eye;
    Vec3f target;
    Vec3f up;
    Vec3f x;
    Vec3f y;
    Vec3f z;
    float aspect;
};

//handle event
void updata_camera_pos(Camera& camera);
void handle_events(Camera& camera);

#endif //RENDERERS_CAMERA_H
