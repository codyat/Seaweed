// Name: Cody Troyer
// Quarter, Year: Spring 2014
// Project 2
//
// This file is to be modified by the student.
// main.cpp
////////////////////////////////////////////////////////////
#include <unistd.h>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <GL/glut.h>

#include "const.h"
#include "color.h"
#include "vector3.h"
#include "particlesystem.h"
#include "object.h"

const int WINDOW_WIDTH = 800;
const int WINDOW_HEIGHT = 800;
const double VIEW_LEFT = 0.0;
const double VIEW_RIGHT = WINDOW_WIDTH;
const double VIEW_BOTTOM = 0.0;
const double VIEW_TOP = WINDOW_HEIGHT;
const double VIEW_FRONT = -800;
const double VIEW_BACK = 800;

int currentTime = 0;
int previousTime = 0;
double tempx, tempy;

std::vector<ParticleSystem*> psystems;
std::vector<Player*> fishies;

Player* me = new Player(Circle2D(Point2D(300, 500), 15, Color(0.0, 0.0, 0.5), true), 0.0, 0.0);

void GLrender();
void GLupdate();
void generate_fishies();
void generate_seaweed();

void key_up(unsigned char key, int x, int y)
{
  switch(key)
  {
    case 'w':
      me->up = false;
      break;
    case 'a':
      me->left = false;
      break;
    case 's':
      me->down = false;
      break;
    case 'd':
      me->right = false;
      break;
    default:
      break;
  }
}

void glKeyboard(unsigned char key, int x, int y)
{
  switch(key)
  {
    case 'w':
      me->up = true;
      break;
    case 'a':
      me->left = true;
      break;
    case 's':
      me->down = true;
      break;
    case 'd':
      me->right = true;
      break;
    case 'p':
      me->vx = 0;
      me->vy = 0;
      tempx = fishies[0]->vx;
      tempy = fishies[0]->vy;
      for(int i = 0; i < fishies.size(); i++)
      {
				fishies[i]->vx = 0;
				fishies[i]->vy = 0;
			}
      break;
    case 'r':
			for(int i = 0; i < fishies.size(); i++)
			{ 
				fishies[i]->vx = tempx;
				fishies[i]->vy = tempy;
			}
			break;
    default:
      break;
  }
}
void collisions()
{
  if(me->c.p.x + me->c.radius > WINDOW_WIDTH)
  {
    me->vx = -me->vx;
    me->c.p.x = WINDOW_WIDTH - me->c.radius;
  }
  else if(me->c.p.x - me->c.radius < 0)
  {
  	me->vx = -me->vx;
    me->c.p.x = 0 + me->c.radius;
  }
  else if(me->c.p.y + me->c.radius > WINDOW_HEIGHT) 
  {
    me->vy = -me->vy;
    me->c.p.y = WINDOW_HEIGHT - me->c.radius;
  }
  else if(me->c.p.y - me->c.radius < 0)
  {
  	me->vy = -me->vy;
    me->c.p.y = 0 + me->c.radius;
  }
  for(int i = 0; i < psystems.size(); i++)
		psystems[i]->collisions(*me);
}

void collisions(int x)
{
  if(fishies[x]->c.p.x + fishies[x]->c.radius > WINDOW_WIDTH)
  {
    fishies[x]->vx = -fishies[x]->vx;
    fishies[x]->c.p.x = WINDOW_WIDTH - fishies[x]->c.radius;
  }
  else if(fishies[x]->c.p.x - fishies[x]->c.radius < 0)
  {
  	fishies[x]->vx = -fishies[x]->vx;
    fishies[x]->c.p.x = 0 + fishies[x]->c.radius;
  }
  else if(fishies[x]->c.p.y + fishies[x]->c.radius > WINDOW_HEIGHT) 
  {
    fishies[x]->vy = -fishies[x]->vy;
    fishies[x]->c.p.y = WINDOW_HEIGHT - fishies[x]->c.radius;
  }
  else if(fishies[x]->c.p.y - fishies[x]->c.radius < 0)
  {
  	fishies[x]->vy = -fishies[x]->vy;
    fishies[x]->c.p.y = 0 + fishies[x]->c.radius;
  }
  for(int i = 0; i < psystems.size(); i++)
		psystems[i]->collisions(*fishies[x]);
}

//Initializes OpenGL attributes
void GLInit(int* argc, char** argv)
{
	glutInit(argc, argv);
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
	glutInitWindowSize(WINDOW_WIDTH, WINDOW_HEIGHT);
	glutCreateWindow("CS 134 - Project 2 - Cody Troyer");
	glutDisplayFunc(GLrender);
	glutIdleFunc(GLupdate);
	glutKeyboardFunc(glKeyboard);
  glutKeyboardUpFunc(key_up);
	glClearColor(0.2f, 0.7f, 1.0f, 0.0f);
	glOrtho(VIEW_LEFT, VIEW_RIGHT, VIEW_BOTTOM, VIEW_TOP, VIEW_FRONT, VIEW_BACK);
	generate_fishies();
	generate_seaweed();
}

int main(int argc, char** argv)
{
	srand(time(NULL));
	GLInit(&argc, argv);
	glutMainLoop();
	return 0;
}

void GLupdate()
{
	const int FRAME_RATE = 25;
	
	double dt = FRAME_RATE / 1000.0;
	updateParticleSystems(psystems, dt);
	cleanupParticleSystems(psystems);

	glutPostRedisplay();

	//sleep is not effective in capturing constant time between frames because sleep
	//doesn't consider the time it takes for context-switching. However, this reduces
	//the cpu-usage. If accurate time frames are desire, use a time accumulator
	currentTime = glutGet(GLUT_ELAPSED_TIME);
	int diffTime = currentTime - previousTime;
	previousTime = currentTime;
	usleep(1000 * std::max(FRAME_RATE - diffTime, 0));
	
	double vx = 0;
	double vy = 0;
	
	if(me->up) vy += .1;
	if(me->down) vy -= .1;
	if(me->right) vx += .1;
	if(me->left) vx -= .1;
	
	me->update_NOW(vx, vy);
	collisions();
	for(int i = 0; i < fishies.size(); i++) 
	{
		collisions(i);
		fishies[i]->update_NOW(0,0);
	}
}

void GLrender()
{
	glClear(GL_COLOR_BUFFER_BIT);
  //print the seaweed behind the fish first
  for(int i = 0; i < psystems.size(); i++) 
		if(!psystems[i]->inFront)
			psystems[i]->render();
	//print the fish next
  for(int i = 0; i < fishies.size(); i++)
    fishies[i]->c.render();
  me->render();
  //now print the seaweed in front of everything
  for(int i = 0; i < psystems.size(); i++) 
		if(psystems[i]->inFront)
			psystems[i]->render();
			
	glFlush();	
	glutSwapBuffers();
}

//generates a random integer between min and max
int randInt(int min, int max)
{
	return (int)(rand() / static_cast<double>(RAND_MAX) * (max - min) + min);
}

//generates a random number of fish
void generate_fishies()
{
	Player *p;
	
	int num = randInt(10,15);
	for(int i = 0; i < num; i++)
	{
		p = new Player(Circle2D(Point2D(randInt(30, 50) * 10, randInt(30, 50) * 10), randInt(3,7), RED, true), 1.5, randInt(1,10)/100.0);
		fishies.push_back(p);
	}
}

// generates a random number of seaweed strands
void generate_seaweed()
{
	ParticleSystemSpringMass* p;
	for(int i = 10; i < WINDOW_WIDTH; i += randInt(40, 80))
	{
		p = new ParticleSystemSpringMass(Vector3(i, 0, 0));
		//randomly decide if the seaweed is infront of the fish, or behind
		p->inFront = rand() % 2;
		psystems.push_back(p);
		if(psystems.size() % randInt(2, 4) == 0) i += randInt(150, 200);
	}
}
