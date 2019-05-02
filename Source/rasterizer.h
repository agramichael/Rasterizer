#ifndef RASTERIZER_H
#define RASTERIZER_H

#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModelH.h"
#include <stdint.h>

using namespace std;
using glm::vec2;
using glm::vec3;
using glm::vec4;
using glm::mat3;
using glm::mat4;
using glm::ivec2;

#define SCREEN_WIDTH 400
#define SCREEN_HEIGHT 400
#define FULLSCREEN_MODE false

SDL_Event event;
vec4 cameraPos( 0, 0, -3.001,1 );
float FOCAL_LENGTH = SCREEN_HEIGHT;
vector<Triangle> triangles;
mat4 R;
mat4 camToWorld;
float yaw = 5 * M_PI / 180;
float total_rot = 0;
float depthBuffer[SCREEN_HEIGHT][SCREEN_WIDTH];
vec3 colorBuffer[SCREEN_HEIGHT][SCREEN_WIDTH];
vec4 lightPos(0,-0.5,-0.7, 1);
vec3 lightPower = 14.f*vec3( 1, 1, 1 );
vec3 indirectLight = 0.5f*vec3( 1, 1, 1 );
vec4 currentNormal;
vec3 currentReflectance;
bool depth_pass = true;
float focal_depth = 2.5;

struct Pixel {
  int x;
  int y;
  float zinv;
  vec4 pos3d;
};

struct Vertex{
  vec4 position;
};

#endif