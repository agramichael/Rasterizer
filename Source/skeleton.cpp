#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModelH.h"
#include <stdint.h>

using namespace std;
using glm::vec2;
using glm::vec3;
using glm::mat3;
using glm::vec4;
using glm::mat4;
using glm::ivec2;


SDL_Event event;

#define SCREEN_WIDTH 400
#define SCREEN_HEIGHT 400
#define FULLSCREEN_MODE false

/* ----------------------------------------------------------------------------*/
/* GLOBALS                                                                   */
vec4 cameraPos( 0, 0, -3.001,1 );
float FOCAL_LENGTH = SCREEN_HEIGHT;
vector<Triangle> triangles;
mat4 R;
mat4 camToWorld;
float yaw = 0.1;

/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

bool Update();
void Draw(screen* screen);
void VertexShader( vec4& v, ivec2& p );
void Interpolate( ivec2 a, ivec2 b, vector<ivec2>& result );
void DrawLineSDL( screen* screen, ivec2 a, ivec2 b, vec3 color );
void DrawPolygonEdges( screen* screen, vector<vec4>& vertices );
mat4 lookAt( vec3 from, vec3 to);

int main( int argc, char* argv[] )
{
  
  screen *screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE );
  LoadTestModel( triangles );
  while ( Update())
    {
      Draw(screen);
      SDL_Renderframe(screen);
    }

  SDL_SaveImage( screen, "screenshot.bmp" );

  KillSDL(screen);
  return 0;
}

/*Place your drawing here*/
void Draw(screen* screen)
{
  /* Clear buffer */
  memset(screen->buffer, 0, screen->height*screen->width*sizeof(uint32_t));
  vec3 color(1,1,1);
 // DrawLineSDL( screen, a,  b,  color );
 for( uint32_t i=0; i<triangles.size(); ++i )
  {
    vector<vec4> vertices(3);
    vertices[0] = triangles[i].v0;
    vertices[1] = triangles[i].v1;
    vertices[2] = triangles[i].v2;
    DrawPolygonEdges(screen, vertices);
  }
}

void update_R(float y) {
  R = mat4(
    cos(y), 0, -sin(y), 0,
    0     , 1, 0      , 0,
    sin(y), 0, cos(y) , 0,
    0     , 0, 0      , 1
  );
}

/*Place updates of parameters here*/
bool Update()
{
  vec4 right( R[0][0], R[0][1], R[0][2], 1 );
  vec4 down( R[1][0], R[1][1], R[1][2], 1 );
  vec4 forward( R[2][0], R[2][1], R[2][2], 1 );

  // static int t = SDL_GetTicks();
  // /* Compute frame time */
  // int t2 = SDL_GetTicks();
  // float dt = float(t2-t);
  // t = t2;

  SDL_Event e;
  while(SDL_PollEvent(&e)) {
    if (e.type == SDL_QUIT) {
      return false;
	  }
  else if (e.type == SDL_KEYDOWN) {
    int key_code = e.key.keysym.sym;
    switch(key_code) {
      case SDLK_UP:
        cameraPos += forward;
	      break;
      case SDLK_DOWN:
        cameraPos -= forward;
    		break;
      case SDLK_LEFT:
        update_R(yaw);
        cameraPos = R * cameraPos;
    		break;
      case SDLK_RIGHT:
	      update_R(-yaw);
        cameraPos = R * cameraPos;
	      break;
      case SDLK_ESCAPE:
	      /* Move camera quit */
	      return false;
      }
	  }  
  }
  return true;
}


void VertexShader( vec4& v, ivec2& p )
{
  mat4 camToWorld = lookAt(vec3(cameraPos), vec3(v));
  v = camToWorld*v;
  p.x = FOCAL_LENGTH * (v.x/v.z) + SCREEN_WIDTH/2;
  p.y = FOCAL_LENGTH * (v.y/v.z) + SCREEN_HEIGHT/2;
}


void Interpolate( ivec2 a, ivec2 b, vector<ivec2>& result )
{
  int N = result.size();
  vec2 step = vec2(b-a) / float(max(N-1,1));
  vec2 current( a );
  for( int i=0; i<N; ++i )
  {
    result[i] = current;
    current += step;
  }
}


void DrawLineSDL( screen* screen, ivec2 a, ivec2 b, vec3 color )
{
  ivec2 delta = glm::abs( a - b );
  int pixels = glm::max( delta.x, delta.y ) + 1;
  vector<ivec2> line( pixels );
  Interpolate( a, b, line );
  for( int i=0; i<line.size(); i++)
  {
    PutPixelSDL( screen, line[i].x, line[i].y, color );
  }
}


void DrawPolygonEdges(screen* screen, vector<vec4>& vertices )
{
  int V = vertices.size();
  // Transform each vertex from 3D world position to 2D image position:
  vector<ivec2> projectedVertices( V );
  for( int i=0; i<V; ++i )
  {
    VertexShader( vertices[i], projectedVertices[i] );
  }
  // Loop over all vertices and draw the edge from it to the next vertex:
  for( int i=0; i<V; ++i )
  {
    int j = (i+1)%V; // The next vertex
    vec3 color( 1, 1, 1 );
    DrawLineSDL( screen, projectedVertices[i], projectedVertices[j], color );
  }
}


mat4 lookAt( vec3 from, vec3 to) {
  vec3 tmp = vec3(0,1,0);
  vec3 forward = normalize(from - to);
  vec3 right = cross(normalize(tmp), forward);
  vec3 up = -cross(forward, right);

  mat4 camToWorld;
  camToWorld[0][0] = right.x;
  camToWorld[0][1] = right.y;
  camToWorld[0][2] = right.z;
  camToWorld[0][3] = 0;

  camToWorld[1][0] = up.x;
  camToWorld[1][1] = up.y;
  camToWorld[1][2] = up.z;
  camToWorld[1][3] = 0;

  camToWorld[2][0] = forward.x;
  camToWorld[2][1] = forward.y;
  camToWorld[2][2] = forward.z;
  camToWorld[2][3] = 0;

  camToWorld[3][0] = from.x;
  camToWorld[3][1] = from.y;
  camToWorld[3][2] = from.z;
  camToWorld[3][3] = 1;

  return camToWorld;
}
