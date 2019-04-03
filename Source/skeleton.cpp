#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>

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
using glm::ivec4;


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
float yaw = 5 * M_PI / 180;
float total_rot = 0;
float depthBuffer[SCREEN_HEIGHT][SCREEN_WIDTH];
vec4 lightPos(0,-0.5,-0.7, 1);
vec3 lightPower = 14.f*vec3( 1, 1, 1 );
vec3 indirectLight = 0.5f*vec3( 1, 1, 1 );
vec4 currentNormal;
vec3 currentReflectance;
struct Pixel {
  int x;
  int y;
  float zinv;
  vec4 pos3d;
};

struct Vertex{
  vec4 position;
};

/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

bool Update();
void Draw(screen* screen);
void VertexShader( Vertex& v, Pixel& p );
void PixelShader( screen* screen, Pixel& p, vec3 color );
void Interpolate( Pixel a, Pixel b, vector<Pixel>& result );
void DrawLineSDL( screen* screen, Pixel a, Pixel b, vec3 color );
void DrawPolygonEdges( screen* screen, vector<vec4>& vertices );
void ComputePolygonRows(  vector<Pixel>& vertexPixels, vector<Pixel>& leftPixels, vector<Pixel>& rightPixels);
void DrawPolygonRows(screen* screen,  vector<Pixel>& leftPixels,  vector<Pixel>& rightPixels, vec3 color);
void DrawPolygon(screen* screen,  vector<Vertex>& vertices, vec3 color );
void update_R(float y);
void load_triangles(std::vector<Triangle>& triangles);
vector<std::string> split(string strToSplit, char delimeter);

int main( int argc, char* argv[] )
{

  screen *screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE );
  load_triangles( triangles );
  for(size_t t = 0; t < triangles.size(); t++){
    cout << triangles[t].v0.x << ", " << triangles[t].v0.y << ", " << triangles[t].v0.z << ", " << endl;
  }
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
  vec3 color;
  for( int y=0; y<SCREEN_HEIGHT; ++y ) {
    for( int x=0; x<SCREEN_WIDTH; ++x ) {
      depthBuffer[y][x] = 0;
    }
  }

 // DrawLineSDL( screen, a,  b,  color );
 for( uint32_t i=0; i<triangles.size(); ++i )
  {
    color = triangles[i].color;
    vector<Vertex> vertices(3);
    vertices[0].position = triangles[i].v0;
    vertices[1].position = triangles[i].v1;
    vertices[2].position = triangles[i].v2;
    currentNormal = triangles[i].normal;
    currentReflectance = 1.f * vec3(1,1,1);
    DrawPolygon(screen, vertices, color);
  }
}

void DrawPolygon(screen* screen, vector<Vertex>& vertices, vec3 color )
{
  int V = vertices.size();
  vector<Pixel> vertexPixels( V );
  for( int i=0; i<V; ++i ) {
    VertexShader( vertices[i], vertexPixels[i] );
  }
  vector<Pixel> leftPixels;
  vector<Pixel> rightPixels;
  ComputePolygonRows( vertexPixels, leftPixels, rightPixels );
  DrawPolygonRows(screen, leftPixels, rightPixels, color );
}

void VertexShader( Vertex& v, Pixel& p )
{
  update_R(total_rot);
  mat4 R_i = glm::inverse(R);
  vec4 current = v.position;
  current = R_i * (current - cameraPos);
  float z = current.z;
  float y = current.y;
  float x = current.x;
  p.zinv = 1 / z;
  p.x = FOCAL_LENGTH * x * p.zinv + SCREEN_WIDTH/2;
  p.y = FOCAL_LENGTH * y * p.zinv + SCREEN_HEIGHT/2;
  p.pos3d = v.position;
}

void ComputePolygonRows( vector<Pixel>& vertexPixels, vector<Pixel>& leftPixels, vector<Pixel>& rightPixels )
{
// 1. Find max and min y-value of the polygon
// and compute the number of rows it occupies.
  int currentMax = -numeric_limits<int>::max();
  int currentMin = +numeric_limits<int>::max();
  for (size_t i = 0; i < vertexPixels.size(); i++){
    if(vertexPixels[i].y > currentMax) currentMax = vertexPixels[i].y;
    if(vertexPixels[i].y < currentMin) currentMin = vertexPixels[i].y;
  }
  int rows = currentMax - currentMin+1;
// 2. Resize leftPixels and rightPixels
// so that they have an element for each row.
  leftPixels.resize(rows);
  rightPixels.resize(rows);
// 3. Initialize the x-coordinates in leftPixels
// to some really large value and the x-coordinates
// in rightPixels to some really small value.
  for( int i=0; i<rows; ++i ) {
    leftPixels[i].x = +numeric_limits<int>::max();
    rightPixels[i].x = -numeric_limits<int>::max();
  }
// 4. Loop through all edges of the polygon and use
// linear interpolation to find the x-coordinate for
// each row it occupies. Update the corresponding
// values in rightPixels and leftPixels.
  for (int i = 0; i < rows; i++) {
    leftPixels[i].y = currentMin + i;
    rightPixels[i].y = currentMin + i;
  }

  int V = vertexPixels.size();
  for( size_t i=0; i<vertexPixels.size(); i++){
    int j = (i+1)%V;
    int rowDiff = abs(vertexPixels[i].y - vertexPixels[j].y);
    vector<Pixel> result(rowDiff + 1);
    Interpolate(vertexPixels[i], vertexPixels[j], result);
    for(size_t k = 0; k < result.size(); k++) {
      int offset =  result[k].y - currentMin;
      if (result[k].x < leftPixels[offset].x) {
        leftPixels[offset].x = result[k].x;
        leftPixels[offset].zinv = result[k].zinv;
        leftPixels[offset].pos3d = result[k].pos3d;
      }
      if (result[k].x > rightPixels[offset].x) {
        rightPixels[offset].x = result[k].x;
        rightPixels[offset].zinv = result[k].zinv;
        rightPixels[offset].pos3d = result[k].pos3d;
      }
    }
  }
}

void DrawPolygonRows(screen* screen,  vector<Pixel>& leftPixels,  vector<Pixel>& rightPixels, vec3 color )
{
  for( size_t row = 0; row < leftPixels.size(); row++ ) {
    DrawLineSDL(screen, leftPixels[row], rightPixels[row], color);
  }
}

void DrawLineSDL( screen* screen, Pixel a, Pixel b, vec3 color )
{
  Pixel delta;
  delta.x = glm::abs( a.x - b.x );
  delta.y = glm::abs( a.y - b.y );
  int pixels = glm::max( delta.x, delta.y ) + 1;
  vector<Pixel> line( pixels );
  Interpolate( a, b, line );
  for( size_t i = 0; i < line.size(); i++){
    PixelShader(screen, line[i], color);
  }
}

void PixelShader( screen* screen, Pixel& p, vec3 color){
  int x = p.x;
  int y = p.y;
  float zinv = p.zinv;
  if( x > 0 && y > 0 && x < SCREEN_WIDTH && y < SCREEN_HEIGHT && zinv >= depthBuffer[y][x] ){
    vec3 r = vec3(lightPos) - vec3(p.pos3d);
    //cout << p.pos3d.x << "," << p.pos3d.y << "," << p.pos3d.z << "," << "\n";
    float radius = glm::length(r);
    r = normalize(r);
    vec3 normal = normalize(vec3(currentNormal));
    float l = max( dot(r, normal), float(0.0) ) / ( 4 * M_PI * pow(radius, 2) );
    vec3 illumination = currentReflectance * (( l * lightPower ) + indirectLight);
    depthBuffer[y][x] = zinv;
    PutPixelSDL( screen, x, y, illumination * color);
  }
}

void update_R(float y) {
  R[0][0] = R[2][2] = cos(y);
  float s = sin(y);
  R[2][0] = s;
  R[0][2] = -s;
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
        total_rot += yaw;
        cameraPos = R * cameraPos;
    		break;
      case SDLK_RIGHT:
	      update_R(-yaw);
        total_rot -= yaw;
        cameraPos = R * cameraPos;
	      break;
      case SDLK_ESCAPE:
	      /* Move camera quit */
	      return false;
      case SDLK_w:
        // move light forwards
        lightPos += 0.2f * forward;
        break;
      case SDLK_s:
        // move light backwards
        lightPos -= 0.2f * forward;
        break;
      case SDLK_a:
        // move light left
        lightPos -= 0.2f * right;
        break;
      case SDLK_d:
        // move light right
        lightPos += 0.2f * right;
        break;
      case SDLK_q:
        // move light up
        lightPos -= 0.2f * down;
        break;
      case SDLK_e:
        // move light down
        lightPos += 0.2f * down;
        break;      
      }
	  }
  }
  return true;
}

void Interpolate( Pixel a, Pixel b, vector<Pixel>& result )
{
  int N = result.size();
  vec4 a_pos = a.pos3d * a.zinv;
  vec4 b_pos = b.pos3d * b.zinv;
  vec3 step;
  step.x = (b.x - a.x) / float(max(N-1,1));
  step.y = (b.y - a.y) / float(max(N-1,1));
  step.z = (b.zinv - a.zinv) / float(max(N-1,1));

  vec4 pos_step = (b_pos - a_pos) / float(max(N-1,1));
  vec3 current;
  vec4 current_position = a_pos;
  current.x = a.x;
  current.y = a.y;
  current.z = a.zinv;
  for( int i=0; i<N; ++i )
  {
    result[i].zinv = current.z;
    result[i].x = round(current.x);
    result[i].y = round(current.y);
    result[i].pos3d = current_position / current.z;
    current_position += pos_step;
    current += step;
  }
}

void DrawPolygonEdges(screen* screen, vector<Vertex>& vertices )
{
  int V = vertices.size();
  // Transform each vertex from 3D world position to 2D image position:
  vector<Pixel> projectedVertices( V );
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



void load_triangles(vector<Triangle>& triangles) {
  string line; //stores
  vector<vec3> vertices;
  vector<ivec4> facets;
  std::ifstream modelFile;
  modelFile.open("Source/cubo.obj", std::ifstream::in);

  if (!modelFile || modelFile.fail()) {
    cout << "Unable to open file" << endl;
    exit(1);   // call system to stop
  }

  // read file line by line create a vec3 for 3D coordinates
  while( std::getline(modelFile, line) ){
    if(line[0] == 'v' && line[1] == ' '){
      stringstream ss(line);
      string item;
      vector<string> splittedStrings;
      int i = 0;
      while(std::getline(ss, item, ' ')){
        if(i != 0){ 
          splittedStrings.push_back(string(item));
        }
        i++;
      }
        vec3 vec(stof(splittedStrings[1]), stof(splittedStrings[2]), stof(splittedStrings[3]));
      vertices.push_back(vec);
    }
    if(line[0] == 'f'){
      stringstream ss(line);
      string item;
      vector<string> splittedStrings;
      int i = 0;
      while(std::getline(ss, item, ' ')){  //1st char in each item = one of facets vertex index.
        if(i != 0) splittedStrings.push_back(string(item));
        i++;
      }
      // split up this thing below with " " and get the first in each
      //splittedStrings[0] == '1/1/1 2/2/2 3/3/3 4/4/4' 
      ivec4 facet;
        for(size_t i = 0; i < splittedStrings.size(); i++){
          vector<string> words = split(splittedStrings[i], ' '); // 1 word = "10/12/11"

          for(size_t j = 0; j < words.size(); j++) { 
            vector<string> items = split(words[j], '/');  //1 item = 10
            if(i != 4)facet[i] = atof(items[0].c_str());
          }
      }
      facets.push_back(facet);
    }
  }
  vec3 green(  0.15f, 0.75f, 0.15f );
  for(size_t i = 0; i < facets.size(); i++){
    vec4 v0 = vec4(vertices[facets[i][0]][0], vertices[facets[i][0]][1], vertices[facets[i][0]][2], 1.0) ;
    vec4 v1 = vec4(vertices[facets[i][1]][0], vertices[facets[i][1]][1], vertices[facets[i][1]][2], 1.0) ;
    vec4 v2 = vec4(vertices[facets[i][2]][0], vertices[facets[i][2]][1], vertices[facets[i][2]][2], 1.0) ;
    vec4 v3 = vec4(vertices[facets[i][3]][0], vertices[facets[i][3]][1], vertices[facets[i][3]][2], 1.0) ;
    triangles.push_back(Triangle( v0, v1, v2, green ));
    triangles.push_back(Triangle( v0, v2, v3, green ));
  }
  modelFile.close();
}


std::vector<std::string> split(std::string strToSplit, char delimeter)
{
    std::stringstream ss(strToSplit);
    std::string item;
    std::vector<std::string> splittedStrings;
    while (std::getline(ss, item, delimeter))
    {
       splittedStrings.push_back(item);
    }
    return splittedStrings;
}