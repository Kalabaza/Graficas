// OpenGL2016P.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "GL/glew.h"
#include "GL/glut.h"
#include <math.h>
#include <algorithm>
#include <windows.h>
#include "Matrix4D\Matrix4D.h"
#include "Mesh\MeshMathSurface.h"
#include "Mesh\ImageBMP.h"
#include <iostream>

GLuint TexFloor;

float time;
bool bFordward, bBackward, bLeft, bRight, bUp, bDown, wireFrame;
bool bTurnLeft, bTurnRight;

bool bVertexSelected;
unsigned long ulVertexSelectedIndex;
MATRIX4D W,V,P,Ac; // Matriz de Mundo, Vista, Proyeccion y Correccion de Aspecto

float mx, my;   // Normalized mouse coordinates

// VECTOR4D Tetraedro[]{ { 0,1,0,1 },{ 1,-1,0,1 },{ -1,-1,0,1 },{ 0,0,1,1 } };
// int TetraedroIndices[]{ 0,1,2,0,2,3,2,1,3,0,3,1 };
// VECTOR4D Colors[]{ { 1,0,0,1 },{ 0,1,0,1 },{ 0,0,1,1 },{ 1,1,1,1 } };

// const float PI = 3.141592654f;
// 
// float SinCos(float x, float y)
// {
//     return 0.3f * cos(3 * x * PI) * sinf(3 * y * PI);
// }

// VECTOR4D SinCosNormal(float x, float y, float z)
// {
//     return Normalize({0.3f * 3 * PI * sinf(3 * x * PI) * sinf(3 * y * PI),
//                      -0.3f * 3 * PI * cosf(3 * x * PI) * cosf(3 * y * PI),
//                       1, 0}); // La componente w es cero
// }

CMesh::VERTEX MyFirstVertexShader(CMesh::VERTEX V)
{
    CMesh::VERTEX out = V;
    VECTOR4D LigthColor{ 1, 1, 1, 1 };
    VECTOR4D LightDirection{ -1, -1, -1, 0 };
    // Iluminacion lambertiana
    float ILambert = max(0.0f, -Dot(V.Normal, LightDirection));
    out.Color = LigthColor * ILambert;
    return out;
}

// float Step(float x, float y)
// {
//     if (x == -1 || y == -1)
//         return 0;
//     float z = -sin(x * PI) * sinf(y * PI);
//     if (z <= 0)
//         return 0;
//     return 0.5;
// }

// float sinc(const float x)
// {
//     if (x == 0)
//         return 1;
//     return sinf(x) / x;
// }

// float WaterDrop(float x, float y)
// {
//     return -10 * sinc(sqrtf(x*x + y*y));
// }
// 
// float Cone(float x, float y)
// {
//     float A = 5.0f, B = 5.0f, C = 10.0f;
//     return -sqrtf(C*C * (((x*x) / (A*A)) + ((y*y) / (B*B))));
// }
// 
// float Plane(float x, float y)
// {
//     float A = 1.0f, B = 1.0f, C = 10.0f, D = 1.0f;
//     return ((D * -1) - (A * x) - (B * y)) / C;
// }
// 
// VECTOR4D Ellipsoid(float u, float v)
// {
//     float A = 1.0f, B = 1.0f, C = 1.0f;
//     VECTOR4D vec{ A * cos(u) * cosf(v), B * cos(u) * sinf(v), C * sinf(u), 1 };
//     return vec;
// }
// 
// VECTOR4D Torus(float u, float v)
// {
//     float R = 1.0f, r = 0.5f;
//     VECTOR4D vec{ (R + r * cos(u)) * cos(v), (R + r * cos(u)) * sinf(v), r * sinf(u), 1 };
//     return vec;
// }
// VECTOR4D Cylinder(float u, float v)
// {
//     float r = 1.0;
//     VECTOR4D vec{ r * cos(u), r * sinf(u), v, 1 };
//     return vec;
// }

CMeshMathSurface g_TheSurface;
// CMeshMathSurface g_Cylinder;
// CMeshMathSurface g_Steps;
// CMeshMathSurface g_WaterDrop;
// CMeshMathSurface g_Cone;
// CMeshMathSurface g_Plane;
// CMeshMathSurface g_Ellipse;
// CMeshMathSurface g_Torus;

void DisplayFunction()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glClearDepth(1.0f);
    int sx, sy;
    sx = glutGet(GLUT_WINDOW_WIDTH);
    sy = glutGet(GLUT_WINDOW_HEIGHT);
    Ac = Scaling(1.0f, static_cast<float>(sx) / static_cast<float>(sy), 1.0f);
    MATRIX4D R = RotationZ(time) * RotationX(time) * RotationY(time);
    P = PerspectiveWidthHeightRH(0.5f, 0.5f, 1.0f, 20.0f);
    //MATRIX4D M = SAspectAc * P * V * R * Translation(0.0f, 0.0f, -1.0f);
    W = Identity();
    MATRIX4D M = Ac * P * V * R * W;
//     MATRIX4D MEllipse = Ac * P * V * R * Scaling(0.8f, 0.8f, 0.8f);
//     MATRIX4D MTorus = Ac * P * V * R * Scaling(0.7f, 0.7f, 0.7f);
//     MATRIX4D MPlane = Ac * P * V * R * Translation(0.0f, 0.0f, 1.1f);
//     MATRIX4D MWater = Ac * P * V * R * Scaling(0.08f, 0.08f, 0.08f);
//     MATRIX4D MCone = Ac * P * V * R * Translation(0.0f, 0.0f, 1.0f) * Scaling(0.7f, 0.7f, 0.7f);
//     MATRIX4D MSteps = Ac * P * V * R * Scaling(0.7f, 0.7f, 0.7f);
//     MATRIX4D MCylinder = Ac * P * V * R * Translation(0.0f, 0.0f, -0.8f) * Scaling(0.7f, 0.7f, 0.7f);
    bool bFill = false;
    multimap<float, CMesh::INTERSECTIONINFO> Faces;
    VECTOR4D WorldRayOrigin, WorldRayDir, ModelRayOrigin, ModelRayDir;
    BuildRayFromPerspective(M, mx, my, WorldRayOrigin, WorldRayDir);
    // Como RayOrigin y RayDir estan en espacio mundial, hay que transformarlos
    // al espacio de modelo donde esta la geometria.
    MATRIX4D InvM;
    Inverse(W, InvM);
    ModelRayOrigin = InvM * WorldRayOrigin;
    ModelRayDir = InvM * WorldRayDir;
    g_TheSurface.RayCast(ModelRayOrigin, ModelRayDir, Faces);
    bFill = Faces.size() != 0 ? true : false;
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, TexFloor);
    glBegin(GL_TRIANGLES);
    g_TheSurface.Draw(M);
    glEnd();
    /*for (int i = 0; i < sizeof(TetraedroIndices) / sizeof(int); i += 3)
    {
        for (int j = 0; j < 3; ++j)
        {
            auto &C = Colors[TetraedroIndices[i + j]];
            glColor4f(C.r, C.g, C.b, C.a);
            VECTOR4D V = M * Tetraedro[TetraedroIndices[i + j]];
            glVertex4f(V.x, V.y, V.z, V.w);
        }
    }*/
//     g_Cylinder.Draw(MCylinder);
//     g_Steps.Draw(MSteps);
//     g_Ellipse.Draw(MEllipse);
//     g_Torus.Draw(MTorus);
//     g_Plane.Draw(MPlane);
//     g_WaterDrop.Draw(MWater);
//     g_Cone.Draw(MCone);
//     glEnd();
    glutSwapBuffers();
}

void IdleFunction()
{
    time += 0.001f;
    MATRIX4D InvV = FastInverse(V);
    VECTOR4D XDir{ InvV.m00, InvV.m10, InvV.m20, 0 };
    VECTOR4D YDir{ InvV.m01, InvV.m11, InvV.m21, 0 };
    VECTOR4D ZDir{ InvV.m02, InvV.m12, InvV.m22, 0 };
    // Parte translativa de la camara
    VECTOR4D EyePos{ InvV.m03, InvV.m13, InvV.m23, 1 };
    VECTOR4D Speed{ 0.005f, 0.005f, 0.005f, 0 };
    MATRIX4D O = InvV;
    O.m03 = 0;
    O.m13 = 0;
    O.m23 = 0;
    if (bFordward)
        EyePos = EyePos - Speed * ZDir;
    if (bBackward)
        EyePos = EyePos + Speed * ZDir;
    if (bTurnLeft)
    {
        MATRIX4D R = RotationAxis(0.005f, YDir);
        O = R * O;
        InvV = O;
    }
    if (bTurnRight)
    {
        MATRIX4D R = RotationAxis(-0.005f, YDir);
        O = R * O;
        InvV = O;
    }
    if (bUp)
        EyePos = EyePos - Speed * YDir;
    if (bDown)
        EyePos = EyePos + Speed * YDir;
    if (bRight)
        EyePos = EyePos - Speed * XDir;
    if (bLeft)
        EyePos = EyePos + Speed * XDir;
    // Solo se cambio la posicion de la camara, por lo que se tiene que actualizar
    InvV.m03 = EyePos.x;
    InvV.m13 = EyePos.y;
    InvV.m23 = EyePos.z;
    V = Orthogonalize(FastInverse(InvV));
    // Produce que se repinte la pantalla
    glutPostRedisplay();
}

void OnKeyDown(unsigned char key, int x, int y)
{
    switch (key)
    {
    // Cuando se presione la tecla Esc se interrumpe la ejecucion del programa
    case 27:
        exit(0);
    case 'a':
        bFordward = true;
        break;
    case 'z':
        bBackward = true;
        break;
    case 'j':
        bTurnLeft = true;
        break;
    case 'l':
        bTurnRight = true;
        break;
    default:
        break;
    }
}

void OnKeyUp(unsigned char key, int x, int y)
{
    switch (key)
    {
    case 'a':
        bFordward = false;
        break;
    case 'z':
        bBackward = false;
        break;
    case 'j':
        bTurnLeft = false;
        break;
    case 'l':
        bTurnRight = false;
        break;
    case 'w':
        if (!wireFrame)
        {
            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
            wireFrame = true;
        }
        else
        {
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
            wireFrame = false;
        }
        break;
    case 's':
        if (bVertexSelected)
        {
            auto &V = g_TheSurface.m_Vertexes[ulVertexSelectedIndex];
            VECTOR4D offset{ V.Normal.x * 0.02f, V.Normal.y * 0.02f, V.Normal.z * 0.02f, 0 };
            V.Position = V.Position + offset;
        }
        break;
    case 'x':
        if (bVertexSelected)
        {
            auto &V = g_TheSurface.m_Vertexes[ulVertexSelectedIndex];
            VECTOR4D offset{ V.Normal.x * 0.02f, V.Normal.y * 0.02f, V.Normal.z * 0.02f, 0 };
            V.Position = V.Position - offset;
        }
        break;
    default:
        break;
    }
}

void OnSpecialDown(int key, int x, int y)
{
    switch (key)
    {
    case GLUT_KEY_DOWN:
        bDown = true;
        break;
    case GLUT_KEY_UP:
        bUp = true;
        break;
    case GLUT_KEY_LEFT:
        bLeft = true;
        break;
    case GLUT_KEY_RIGHT:
        bRight = true;
        break;
    default:
        break;
    }
}

void OnSpecialUp(int key, int x, int y)
{
    switch (key)
    {
    case GLUT_KEY_DOWN:
        bDown = false;
        break;
    case GLUT_KEY_UP:
        bUp = false;
        break;
    case GLUT_KEY_LEFT:
        bLeft = false;
        break;
    case GLUT_KEY_RIGHT:
        bRight = false;
        break;
    default:
        break;
    }
}

void OnMouse(int button, int state, int x, int y)
{
    // button:=(GLUT_LEFT_BUTTON|GLUT_MIDDLE_BUTTON|GLUT_RIGHT_BUTTON)
    // state:=(GLUT_DOWN|GLUT_UP)
    mx = -1 + 2 * (float)x / glutGet(GLUT_WINDOW_WIDTH);
    my = 1 - 2 * (float)y / glutGet(GLUT_WINDOW_HEIGHT);
    if (state == GLUT_DOWN)
    {
        VECTOR4D RayOrigin, RayDirection;
        BuildRayFromPerspective(Ac*P*V*W, mx, my, RayOrigin, RayDirection);
        multimap<float, unsigned long> Vertices;
        if (g_TheSurface.RayCast(RayOrigin, RayDirection, Vertices, 0.05f))
        {
            bVertexSelected = true;
            ulVertexSelectedIndex = Vertices.begin()->second;
        }
    }
    if (state == GLUT_UP)
    {
        bVertexSelected = false;
    }
}

void OnMotion(int x, int y)
{
    mx = -1 + 2 * (float)x / glutGet(GLUT_WINDOW_WIDTH);
    my = 1 - 2 * (float)y / glutGet(GLUT_WINDOW_HEIGHT);
}

int main(int argc, char* argv[])
{
//     VECTOR4D a{ 0, 1, 0, 1 };
//     VECTOR4D b{ 1, -1, 0, 1 };
//     VECTOR4D c{ -1, -1, 0, 1 };
// 
//     VECTOR4D p1{ 0.0f, -1.1f, 0.0f, 1.0f };
//     VECTOR4D p2{ 0.0f, -0.9f, 0.0f, 1.0f };
// 
//     cout << (PtInTriangle(a, b, c, p1) == true ? "Dentro" : "Fuera") << endl;   // Fuera
//     cout << (PtInTriangle(a, b, c, p2) == true ? "Dentro" : "Fuera") << endl;   // Dentro
//     cout << (PtInTriangle(a, b, c, a) == true ? "Dentro" : "Fuera") << endl;    // Fuera

//     float w0, w1, w2;
// 
//     cout << (PtInTriangleBarycentric(a, b, c, p1, w0, w1, w2) == true ? "Dentro" : "Fuera") << endl;   // Fuera
//     cout << (PtInTriangleBarycentric(a, b, c, p2, w0, w1, w2) == true ? "Dentro" : "Fuera") << endl;   // Dentro
//     cout << (PtInTriangleBarycentric(a, b, c, a, w0, w1, w2) == true ? "Dentro" : "Fuera") << endl;    // Dentro

    // Inicializar la vista una sola vez al iniciar la ejecucion del programa
    VECTOR4D Target{ 0, 0, 0, 1 };
    VECTOR4D Eye{ 2, 2, 2, 1 };
    VECTOR4D Up{ 0, 0, 1, 0 };
    V = View(Eye, Target, Up);
    
//     unsigned long steps = 25;
    VECTOR4D A{ 1, 0, 0, 1 }, B{ 0, 1, 0, 1 }, C{ 0, 0, 1, 1 }, D{ 1, 1, 1, 1 };

    // Generacion de una malla por medio de una funcion analitica
    g_TheSurface.BuildCube();
//     g_TheSurface.SetColor(A, B, C, D);
//     g_TheSurface.BuildAnalyticSurface(steps, steps, -1, -1, 2.0f / (steps - 1), 2.0f / (steps - 1), SinCos, SinCosNormal);
    // Aplicacion de colores a la malla generada
    g_TheSurface.VertexShade(MyFirstVertexShader);
//     g_TheSurface.BuildTextureCoords(0, 0, 3.0f / (steps - 1), 3.0f / (steps - 1));
//     g_Steps.BuildAnalyticSurface(steps - 1 , steps - 1, -1, -1, 2.0f / (steps - 1), 2.0f / (steps - 1), Step);
//     g_Steps.SetColor(A, B, C, D);
//     g_Cone.BuildAnalyticSurface(steps, steps, -1, -1, 2.0f / (steps - 1), 2.0f / (steps - 1), Cone);
//     g_Cone.SetColor(A, B, C, D);
//     g_WaterDrop.BuildAnalyticSurface(100, 100, -10, -10, 2.0f / (10 - 1), 2.0f / (10 - 1), WaterDrop);
//     g_WaterDrop.SetColor(A, B, C, D);
//     g_Plane.BuildAnalyticSurface(steps, steps, -1, -1, 1.0f / (steps - 1), 1.0f / (steps - 1), Plane);
//     g_Plane.SetColor(A, B, C, D);
//     // Generacion de una malla para una elipse/esfera
//     g_Ellipse.BuildParametricSurface(steps, steps, -PI / 2, -PI, (2 * PI) / (steps - 1), PI / (steps -1), Ellipsoid);
//     g_Ellipse.SetColor(A, B, C, D);
//     // Generacion de una malla para un toroide
//     g_Torus.BuildParametricSurface(steps, steps, 0, 0, (2 * PI) / (steps - 1), (2 * PI) / (steps - 1), Torus);
//     g_Torus.SetColor(A, B, C, D);
//     g_Cylinder.BuildParametricSurface(steps, steps, 0, -1, 2.0f * PI / (steps - 1), 2.0f / (steps - 1), Cylinder);
//     g_Cylinder.SetColor(A, B, C, D);
    // glutinit inicializa la biblioteca de glut
    glutInit(&argc, argv);
    // propiedades de la ventana
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);
    // Tamaño de la ventana
    glutInitWindowSize(1024, 768);
    // Posicion inicial en donde aparece la ventana
    glutInitWindowPosition(100, 100);
    // Crear la ventana de OpenGL
    glutCreateWindow("Mi primer programa en OpenGL");

    glClearColor(0.5f, 0.5f, 0.5f, 0.0f);

    // Despues de crear la ventana se carga la textura en OpenGL
    //CImageBMP *pImage = CImageBMP::CreateBitmapFromFile("..\\Resources\\cube.bmp", NULL);
    //CImageBMP *pImage = CImageBMP::CreateBitmapFromFile("..\\Resources\\dice.bmp", NULL);
    CImageBMP *pImage = CImageBMP::CreateBitmapFromFile("..\\Resources\\minecraft.bmp", NULL);
    TexFloor = pImage->CreateTexure();
    CImageBMP::DestroyBitmap(pImage);

    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);
    glDepthRange(0, 1);
    // Callback que se llamara cada que se ejecuta una iteracion en el ciclo principal
    glutDisplayFunc(DisplayFunction);
    // Callback para la funcion que se ejecutara cuando este en reposo (idle)
    glutIdleFunc(IdleFunction);
    glutKeyboardFunc(OnKeyDown);
    glutKeyboardUpFunc(OnKeyUp);
    glutSpecialFunc(OnSpecialDown);
    glutSpecialUpFunc(OnSpecialUp);
    glutMouseFunc(OnMouse);
    glutPassiveMotionFunc(OnMotion);
    // Ciclo principal
    glutMainLoop();

    return 0;
}
