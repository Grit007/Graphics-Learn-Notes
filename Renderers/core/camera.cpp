#include "camera.h"
#include "../platform/win32.h"

Camera::Camera(Vec3f e, Vec3f t, Vec3f up, float aspect):
        eye(e),target(t),up(up),aspect(aspect)
{}

Camera::~Camera()
{}

void updata_camera_pos(Camera &camera)
{
    Vec3f from_target = camera.eye - camera.target;
    float radius = from_target.norm();

    float phi     = (float)atan2(from_target[0], from_target[2]); // azimuth angle(方位角), angle between from_target and z-axis，[-pi, pi]
    float theta   = (float)acos(from_target[1] / radius);		  // zenith angle(天顶角), angle between from_target and y-axis, [0, pi]
    float x_delta = window->mouse_info.orbit_delta[0] / window->width;
    float y_delta = window->mouse_info.orbit_delta[1] / window->height;

    // 滑动鼠标事件
    radius *= (float)pow(0.95, window->mouse_info.wheel_delta);

    float factor = 1.5 * PI;
    // 鼠标左键
    phi	  += x_delta * factor;
    theta += y_delta * factor;
    if (theta > PI) theta = PI - EPSILON * 100;
    if (theta < 0)  theta = EPSILON * 100;

    camera.eye[0] = camera.target[0] + radius * sin(phi) * sin(theta);
    camera.eye[1] = camera.target[1] + radius * cos(theta);
    camera.eye[2] = camera.target[2] + radius * sin(theta) * cos(phi);

    // 鼠标右键
    factor  = radius * (float)tan(60.0 / 360 * PI) * 2.2;
    x_delta = window->mouse_info.fv_delta[0] / window->width;
    y_delta = window->mouse_info.fv_delta[1] / window->height;
    Vec3f left = camera.x * x_delta * factor;
    Vec3f up   = camera.y * y_delta * factor;

    camera.eye += (left - up);
    camera.target += (left - up);
}

void handle_mouse_events(Camera& camera)
{
    // 注意 鼠标事件是当鼠标按下并移动的时候才会work
    if (window->buttons[0])
    {
        Vec2f cur_pos = get_mouse_pos();
        window->mouse_info.orbit_delta = window->mouse_info.orbit_pos - cur_pos;
        window->mouse_info.orbit_pos = cur_pos;
    }

    if (window->buttons[1])
    {
        Vec2f cur_pos = get_mouse_pos();
        window->mouse_info.fv_delta = window->mouse_info.fv_pos - cur_pos;
        window->mouse_info.fv_pos = cur_pos;
    }

    updata_camera_pos(camera);
}
extern bool blinn;
void handle_key_events(Camera& camera)
{
    float distance = (camera.target - camera.eye).norm();
    // 是否开启blinn着色
    if (window->keys['B'])
    {
        blinn = !blinn;
    }
    // 放大
    if (window->keys['W'])
    {
        camera.eye += camera.z * distance * (-10.0 / window->width);
    }
    // 缩小
    if (window->keys['S'])
    {
        camera.eye += camera.z * distance * (10.0 / window->width);
    }
    // 向上移动
    if (window->keys[VK_UP] || window->keys['Q'])
    {
        camera.eye -= camera.y * 0.05f;
        camera.target -= camera.y * 0.05f;
    }
    // 向下移动
    if (window->keys[VK_DOWN] || window->keys['E'])
    {
        camera.eye -= camera.y * (-0.05f);
        camera.target -= camera.y * (-0.05f);
    }
    // 向左移动
    if (window->keys[VK_LEFT] || window->keys['A'])
    {
        camera.eye -= camera.x * (-0.05f);
        camera.target -= camera.x * (-0.05f);
    }
    // 向右移动
    if (window->keys[VK_RIGHT] || window->keys['D'])
    {
        camera.eye -= camera.x * 0.05f;
        camera.target -= camera.x * 0.05f;
    }
    // ESC关闭窗口
    if (window->keys[VK_ESCAPE])
    {
        window->is_close = 1;
    }
}

void handle_events(Camera& camera)
{
    //calculate camera axis
    camera.z = (camera.eye - camera.target).normalize();
    camera.x = (cross(camera.up, camera.z)).normalize();
    camera.y = (cross(camera.z, camera.x)).normalize();

    //mouse and keyboard events
    handle_mouse_events(camera);
    handle_key_events(camera);
}
