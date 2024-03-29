#include "rasterizer.h"

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
void FilterCreation(float GKernel[3][3], float sigma);

int main() {
  screen *screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE );
  LoadTestModel( triangles );
  update_R(0);

  while ( Update() ) {
      Draw(screen);
      SDL_Renderframe(screen);
    }

  SDL_SaveImage(screen, "screenshot.bmp");
  KillSDL(screen);
  return 0;
}

// Render a frame
void Draw(screen* screen) {
  // Clear buffers
  memset(screen->buffer, 0, screen->height*screen->width*sizeof(uint32_t));
  memset(depthBuffer, 0, screen->height*screen->width*sizeof(float));
  memset(colorBuffer, 0, screen->height*screen->width*sizeof(vec3));

  vec3 color;

  // Depth pass - do one pass without drawing, only for collecting depth information
  depth_pass = true;
  #pragma omp parallel for schedule(static,10)
  for( uint32_t i=0; i<triangles.size(); ++i ) {
    color = triangles[i].color;
    vector<Vertex> vertices(3);
    vertices[0].position = triangles[i].v0;
    vertices[1].position = triangles[i].v1;
    vertices[2].position = triangles[i].v2;
    currentNormal = triangles[i].normal;
    currentReflectance = 1.f * vec3(1,1,1);
    DrawPolygon(screen, vertices, color);
  }

  // Normal pass - for every triangle, store its vertices, normal and color: use this to draw them.
  depth_pass = false;
  #pragma omp parallel for schedule(static,10)
  for( uint32_t i=0; i<triangles.size(); ++i ) {
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

// Draw a polygon
void DrawPolygon(screen* screen, vector<Vertex>& vertices, vec3 color )
{
  int V = vertices.size();
  vector<Pixel> vertexPixels( V );
  // Project each vertex into the image space
  for( int i=0; i<V; ++i ) {
    VertexShader( vertices[i], vertexPixels[i] );
  }
  vector<Pixel> leftPixels;
  vector<Pixel> rightPixels;
  // Compute beginning and ending of each row
  ComputePolygonRows( vertexPixels, leftPixels, rightPixels );
  // Draw rows
  DrawPolygonRows(screen, leftPixels, rightPixels, color );
}

// Project a vertex into image space: save the inverse of its depth and its original 3d position
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

// Compute beginning and ending of each row in the polygon
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

// Draw the rows of a polygon
void DrawPolygonRows(screen* screen,  vector<Pixel>& leftPixels,  vector<Pixel>& rightPixels, vec3 color )
{
  for( size_t row = 0; row < leftPixels.size(); row++ ) {
    DrawLineSDL(screen, leftPixels[row], rightPixels[row], color);
  }
}

// Draw a line from a to b (interpolating in between)
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

// Draw one pixel: calculating illumination and blur factors
void PixelShader( screen* screen, Pixel& p, vec3 color){
  int x = p.x;
  int y = p.y;
  float zinv = p.zinv;
  if( x > 0 && y > 0 && x < SCREEN_WIDTH && y < SCREEN_HEIGHT && zinv >= depthBuffer[y][x]){
    vec3 r = vec3(lightPos) - vec3(p.pos3d);
    //cout << p.pos3d.x << "," << p.pos3d.y << "," << p.pos3d.z << "," << "\n";
    float radius = glm::length(r);
    r = normalize(r);
    vec3 normal = normalize(vec3(currentNormal));
    float l = max( dot(r, normal), float(0.0) ) / ( 4 * M_PI * pow(radius, 2) );
    vec3 illumination = currentReflectance * (( l * lightPower ) + indirectLight);
    if (depth_pass) {
      depthBuffer[y][x] = zinv;
      colorBuffer[y][x] = color;
    }
    else{
      float current_depth = 1/depthBuffer[y][x];
      float kernel[3][3];
      float diff = abs(focal_depth - current_depth);
      FilterCreation(kernel, diff);

      vec3 blur = vec3(0,0,0);

      for (int i = -1; i < 2; i++) {
        for (int j = -1; j < 2; j++) {
          blur += kernel[i+1][j+1] * colorBuffer[y + i][x + j];
        }
      }

      // cout << "Blur vec3: " << blur.x << ", " << blur.y << ", " << blur.z << "\n";
      PutPixelSDL( screen, x, y, illumination * blur);
    }
  }
}

void FilterCreation(float GKernel[][3], float sigma) { 
    float r, s = 2.0 * sigma * sigma; 
    float sum = 0.0; 
  
    // Generating 3x3 kernel 
    for (int x = -1; x <= 1; x++) { 
        for (int y = -1; y <= 1; y++) { 
            r = sqrt(x * x + y * y); 
            GKernel[x + 1][y + 1] = (exp(-(r * r) / s)) / (M_PI * s); 
            sum += GKernel[x + 1][y + 1]; 
        } 
    } 
  
    // Normalising the kernel 
    for (int i = 0; i < 3; ++i) 
        for (int j = 0; j < 3; ++j) 
            GKernel[i][j] /= sum; 
} 

// Update the rotation matrix
void update_R(float y) {
  R[0][0] = R[2][2] = cos(y);
  float s = sin(y);
  R[2][0] = s;
  R[0][2] = -s;
}

// Camera and light movement with rotation matrix updates
bool Update()
{
  vec4 right( R[0][0], R[0][1], R[0][2], 1 );
  vec4 down( R[1][0], R[1][1], R[1][2], 1 );
  vec4 forward( R[2][0], R[2][1], R[2][2], 1 );

  SDL_Event e;
  while(SDL_PollEvent(&e)) {
    if (e.type == SDL_QUIT) {
      return false;
	  }
  else if (e.type == SDL_KEYDOWN) {
    int key_code = e.key.keysym.sym;
    switch(key_code) {
      case SDLK_UP:
        cameraPos += 0.2f * forward;
	      break;
      case SDLK_DOWN:
        cameraPos -= 0.2f * forward;
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

// Interpolate from pixel a to pixel b
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

// Unused function which draws the edges of a polygon
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
