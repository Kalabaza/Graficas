// OpenGL2016P.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "GL/glew.h"
#include "GL/glut.h"
#include <math.h>
#include <windows.h>
#include "Matrix4D\Matrix4D.h"

float time;

VECTOR4D Tetraedro[]{ { 0,1,0,1 },{ 1,-1,0,1 },{ -1,-1,0,1 },{ 0,0,1,1 } };
int TetraedroIndices[]{ 0,1,2,0,2,3,2,1,3,0,3,1 };
VECTOR4D Colors[]{ { 1,0,0,1 },{ 0,1,0,1 },{ 0,0,1,1 },{ 1,1,0,1 } };

void DisplayFunction()
{
    time += 0.01f;
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glClearDepth(1.0f);
    int sx, sy;
    sx = glutGet(GLUT_WINDOW_WIDTH);
    sy = glutGet(GLUT_WINDOW_HEIGHT);
    MATRIX4D SAspect = Scaling(1, static_cast<float>(sx) / static_cast<float>(sy), 1);
    MATRIX4D T = Translation(0, 0, 0);
    VECTOR4D Target{ 0, 0, 0, 1 };
    VECTOR4D Eye{ 0.1f, 0.1f, 0.1f, 1 };
    VECTOR4D Up{ 0, 0, 1, 0 };
    MATRIX4D V = View(Eye, Target, Up);
    MATRIX4D R = SAspect * V * RotationZ(time) * Scaling(0.5, 0.5, 0.5);
    glBegin(GL_TRIANGLES);
    for (int i = 0; i < sizeof(TetraedroIndices) / sizeof(int); i += 3)
    {
        for (int j = 0; j < 3; ++j)
        {
            auto &C = Colors[TetraedroIndices[i + j]];
            glColor4f(C.r, C.g, C.b, C.a);
            VECTOR4D V = R*Tetraedro[TetraedroIndices[i + j]];
            glVertex4f(V.x, V.y, V.z, V.w);
        }
    }
    glEnd();
    glutSwapBuffers();
}

int main(int argc, char* argv[])
{
    // glutinit inicializa la biblioteca de glut
    glutInit(&argc, argv);
    // propiedades de la ventana
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);
    // tamaño de la ventana
    glutInitWindowSize(400, 400);
    // posicion inicial en donde aparece la ventana
    glutInitWindowPosition(100, 100);
    //glEnable(GL_DEPTH_TEST);
    // crear la ventana de OpenGL
    glutCreateWindow("Mi primer programa en OpenGL");
    glEnable(GL_DEPTH_TEST);
    glDepthRange(0, 1);
    // callback que se llamara cada que se ejecuta una iteracion en el ciclo principal
    glutDisplayFunc(DisplayFunction);
    // ciclo principal
    glutMainLoop();

    return 0;
}
