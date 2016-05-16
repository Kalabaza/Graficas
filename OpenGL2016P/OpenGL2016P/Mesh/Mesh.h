#pragma once

#include "..\Matrix4D\Matrix4D.h"
#include <vector>
#include <map>

using namespace std;

class CMesh
{
public:
    struct VERTEX
    {
        VECTOR4D Position;
        VECTOR4D Normal;
        VECTOR4D Color;
        VECTOR4D TexCoord;
    };
    struct INTERSECTIONINFO
    {
        int Face;
        VECTOR4D LocalPosition;
    };
    // Arreglo de vertices
    vector<VERTEX> m_Vertexes;
    // Buffer de indices
    vector<unsigned int> m_Indexes;
    
    CMesh();
    ~CMesh();

    // Metodo para dibujar la maya en pantalla
    void Draw(const MATRIX4D &M);
    void Draw(const MATRIX4D &M, int StartFace, int Faces);

    bool RayCast(VECTOR4D &RayOrigin, VECTOR4D &RayDir, multimap<float, INTERSECTIONINFO> &Faces);
    bool RayCast(VECTOR4D &RayOrigin, VECTOR4D &RayDir, multimap<float, unsigned long> &Vertices, float radius);
    void VertexShade(VERTEX(*pVS)(VERTEX V));

    void BuildCube();
};

bool RaySphereIntersect(VECTOR4D &RayOrigin, VECTOR4D &RayDirection, VECTOR4D &SphereCenter, float r);