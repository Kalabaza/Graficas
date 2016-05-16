#include "../stdafx.h"
#include "Mesh.h"
#include "GL/glew.h"

CMesh::CMesh()
{
}

CMesh::~CMesh()
{
}

void CMesh::Draw(const MATRIX4D &M)
{
    for (unsigned int i = 0; i < m_Indexes.size(); i += 3)
    {
        for (int j = 0; j < 3; ++j)
        {
            auto &C = m_Vertexes[m_Indexes[i + j]].Color;
            glColor4f(C.r, C.g, C.b, C.a);
            auto &T = m_Vertexes[m_Indexes[i + j]].TexCoord;
            glTexCoord2f(T.x, T.y);
            VECTOR4D V = M * m_Vertexes[m_Indexes[i + j]].Position;
            glVertex4f(V.x, V.y, V.z, V.w);
        }
    }
}

void CMesh::Draw(const MATRIX4D &M, int StartFace, int Faces)
{
    int nTotalFaces = m_Indexes.size() / 3;
    if (Faces == -1)
        Faces = nTotalFaces / 3;
    if (StartFace + Faces > nTotalFaces)
        Faces = nTotalFaces - StartFace;
    for (unsigned int i = StartFace * 3; Faces--; i += 3)
    {
        for (int j = 0; j < 3; ++j)
        {
            auto &C = m_Vertexes[m_Indexes[i + j]].Color;
            glColor4f(C.r, C.g, C.b, C.a);
            auto &T = m_Vertexes[m_Indexes[i + j]].TexCoord;
            glTexCoord2f(T.x, T.y);
            VECTOR4D V = M * m_Vertexes[m_Indexes[i + j]].Position;
            glVertex4f(V.x, V.y, V.z, V.w);
        }
    }
}

bool CMesh::RayCast(VECTOR4D &RayOrigin, VECTOR4D &RayDir, multimap<float, INTERSECTIONINFO> &Faces)
{
    unsigned long nFaces = m_Indexes.size() / 3;
    unsigned long nBaseIndex = 0;
    unsigned long nIntersectedFaces = 0;
    for (unsigned long iFace = 0; iFace < nFaces; ++iFace, nBaseIndex+=3)
    {
        VECTOR4D &V0 = m_Vertexes[m_Indexes[nBaseIndex]].Position;
        VECTOR4D &V1 = m_Vertexes[m_Indexes[nBaseIndex + 1]].Position;
        VECTOR4D &V2 = m_Vertexes[m_Indexes[nBaseIndex + 2]].Position;
        VECTOR4D Intersection;
        if (RayCastOnTriangle(V0, V1, V2, RayOrigin, RayDir, Intersection))
        {
            // Distancia entre el origen y esa interseccion
            float dist = Magnity(Intersection - RayOrigin);
            INTERSECTIONINFO II;
            II.Face = iFace;
            II.LocalPosition = Intersection;
            Faces.insert(make_pair(dist, II));
            ++nIntersectedFaces;
        }
    }
    return nIntersectedFaces != 0;
}

bool RaySphereIntersect(VECTOR4D &RayOrigin, VECTOR4D &RayDirection, VECTOR4D &SphereCenter, float r)
{
    // V^2 = Dot(V, V)
    VECTOR4D RO = RayOrigin - SphereCenter;
    float a = Dot(RayDirection, RayDirection);
    float b = 2 * Dot(RayDirection, RO);
    float c = Dot(RO, RO) - r * r;
    float disc = b * b - 4 * a * c;
    if (disc < 0)
        return false;
    return true;
}

bool CMesh::RayCast(VECTOR4D &RayOrigin, VECTOR4D &RayDir, multimap<float, unsigned long> &Vertices, float radius)
{
    // Ecuacion de la esfera con centro en el origen
    // x^2 + y^2 + z^2 = r^2
    for (unsigned long index = 0; index < m_Vertexes.size(); ++index)
    {
        if (RaySphereIntersect(RayOrigin, RayDir, m_Vertexes[index].Position, radius))
        {
            Vertices.insert(make_pair(Magnity(m_Vertexes[index].Position - RayOrigin), index));
        }
    }
    return Vertices.size() != 0;
}

void CMesh::VertexShade(VERTEX(*pVS)(VERTEX V))
{
    for (size_t i = 0; i < m_Vertexes.size(); ++i)
        m_Vertexes[i] = pVS(m_Vertexes[i]);
}

void CMesh::BuildCube()
{
    VECTOR4D P1{ 0.5, -0.5, -0.5, 1.0 };
    VECTOR4D P2{ 0.5, -0.5,  0.5, 1.0 };
    VECTOR4D P3{ -0.5, -0.5,  0.5, 1.0 };
    VECTOR4D P4{ -0.5, -0.5, -0.5, 1.0 };
    VECTOR4D P5{ 0.5,  0.5, -0.5, 1.0 };
    VECTOR4D P6{ 0.5,  0.5,  0.5, 1.0 };
    VECTOR4D P7{ -0.5,  0.5,  0.5, 1.0 };
    VECTOR4D P8{ -0.5,  0.5, -0.5, 1.0 };

    VECTOR4D T1{ 0.748573f, 0.750412f, 0.0, 0.0 };
    VECTOR4D T2{ 0.749279f, 0.501284f, 0.0, 0.0 };
    VECTOR4D T3{ 0.999110f, 0.501077f, 0.0, 0.0 };
    VECTOR4D T4{ 0.999455f, 0.750380f, 0.0, 0.0 };
    VECTOR4D T5{ 0.250471f, 0.500702f, 0.0, 0.0 };
    VECTOR4D T6{ 0.249682f, 0.749677f, 0.0, 0.0 };
    VECTOR4D T7{ 0.001085f, 0.750380f, 0.0, 0.0 };
    VECTOR4D T8{ 0.001517f, 0.499994f, 0.0, 0.0 };
    VECTOR4D T9{ 0.499422f, 0.500239f, 0.0, 0.0 };
    VECTOR4D T10{ 0.500149f, 0.750166f, 0.0, 0.0 };
    VECTOR4D T11{ 0.748355f, 0.998230f, 0.0, 0.0 };
    VECTOR4D T12{ 0.500193f, 0.998728f, 0.0, 0.0 };
    VECTOR4D T13{ 0.498993f, 0.250415f, 0.0, 0.0 };
    VECTOR4D T14{ 0.748953f, 0.250920f, 0.0, 0.0 };

    // Normal de 1 para suprimir el efecto de la luz
    VECTOR4D N1{ 1.0f, 1.0f, 1.0f, 1.0 };
    VECTOR4D N2{ 1.0f, 1.0f, 1.0f, 1.0 };
    VECTOR4D N3{ 1.0f, 1.0f, 1.0f, 1.0 };
    VECTOR4D N4{ 1.0f, 1.0f, 1.0f, 1.0 };
    VECTOR4D N5{ 1.0f, 1.0f, 1.0f, 1.0 };
    VECTOR4D N6{ 1.0f, 1.0f, 1.0f, 1.0 };
    VECTOR4D N7{ 1.0f, 1.0f, 1.0f, 1.0 };
    VECTOR4D N8{ 1.0f, 1.0f, 1.0f, 1.0 };

    //     VECTOR4D N1{  0.000000f,  0.000000f, -1.000000f, 1.0};
    //     VECTOR4D N2{ -1.000000f, -0.000000f, -0.000000f, 1.0 };
    //     VECTOR4D N3{ -0.000000f, -0.000000f,  1.000000f, 1.0 };
    //     VECTOR4D N4{ -0.000001f,  0.000000f,  1.000000f, 1.0 };
    //     VECTOR4D N5{  1.000000f, -0.000000f,  0.000000f, 1.0 };
    //     VECTOR4D N6{  1.000000f,  0.000000f,  0.000001f, 1.0 };
    //     VECTOR4D N7{  0.000000f,  1.000000f, -0.000000f, 1.0 };
    //     VECTOR4D N8{ -0.000000f, -1.000000f,  0.000000f, 1.0 };

    m_Vertexes.resize(36);
    m_Vertexes[0].Position = P5; m_Vertexes[1].Position = P1; m_Vertexes[2].Position = P4;
    m_Vertexes[3].Position = P5; m_Vertexes[4].Position = P4; m_Vertexes[5].Position = P8;
    m_Vertexes[6].Position = P3; m_Vertexes[7].Position = P7; m_Vertexes[8].Position = P8;
    m_Vertexes[9].Position = P3; m_Vertexes[10].Position = P8; m_Vertexes[11].Position = P4;
    m_Vertexes[12].Position = P2; m_Vertexes[13].Position = P6; m_Vertexes[14].Position = P3;
    m_Vertexes[15].Position = P6; m_Vertexes[16].Position = P7; m_Vertexes[17].Position = P3;
    m_Vertexes[18].Position = P1; m_Vertexes[19].Position = P5; m_Vertexes[20].Position = P2;
    m_Vertexes[21].Position = P5; m_Vertexes[22].Position = P6; m_Vertexes[23].Position = P2;
    m_Vertexes[24].Position = P5; m_Vertexes[25].Position = P8; m_Vertexes[26].Position = P6;
    m_Vertexes[27].Position = P8; m_Vertexes[28].Position = P7; m_Vertexes[29].Position = P6;
    m_Vertexes[30].Position = P1; m_Vertexes[31].Position = P2; m_Vertexes[32].Position = P3;
    m_Vertexes[33].Position = P1; m_Vertexes[34].Position = P3; m_Vertexes[35].Position = P4;

    m_Vertexes[0].TexCoord = T1; m_Vertexes[1].TexCoord = T2; m_Vertexes[2].TexCoord = T3;
    m_Vertexes[3].TexCoord = T1; m_Vertexes[4].TexCoord = T3; m_Vertexes[5].TexCoord = T4;
    m_Vertexes[6].TexCoord = T5; m_Vertexes[7].TexCoord = T6; m_Vertexes[8].TexCoord = T7;
    m_Vertexes[9].TexCoord = T5; m_Vertexes[10].TexCoord = T7; m_Vertexes[11].TexCoord = T8;
    m_Vertexes[12].TexCoord = T9; m_Vertexes[13].TexCoord = T10; m_Vertexes[14].TexCoord = T5;
    m_Vertexes[15].TexCoord = T10; m_Vertexes[16].TexCoord = T6; m_Vertexes[17].TexCoord = T5;
    m_Vertexes[18].TexCoord = T2; m_Vertexes[19].TexCoord = T1; m_Vertexes[20].TexCoord = T9;
    m_Vertexes[21].TexCoord = T1; m_Vertexes[22].TexCoord = T10; m_Vertexes[23].TexCoord = T9;
    m_Vertexes[24].TexCoord = T1; m_Vertexes[25].TexCoord = T11; m_Vertexes[26].TexCoord = T10;
    m_Vertexes[27].TexCoord = T11; m_Vertexes[28].TexCoord = T12; m_Vertexes[29].TexCoord = T10;
    m_Vertexes[30].TexCoord = T2; m_Vertexes[31].TexCoord = T9; m_Vertexes[32].TexCoord = T13;
    m_Vertexes[33].TexCoord = T2; m_Vertexes[34].TexCoord = T13; m_Vertexes[35].TexCoord = T14;

    m_Vertexes[0].Normal = N1; m_Vertexes[1].Normal = N1; m_Vertexes[2].Normal = N1;
    m_Vertexes[3].Normal = N1; m_Vertexes[4].Normal = N1; m_Vertexes[5].Normal = N1;
    m_Vertexes[6].Normal = N2; m_Vertexes[7].Normal = N2; m_Vertexes[8].Normal = N2;
    m_Vertexes[9].Normal = N2; m_Vertexes[10].Normal = N2; m_Vertexes[11].Normal = N2;
    m_Vertexes[12].Normal = N3; m_Vertexes[13].Normal = N3; m_Vertexes[14].Normal = N3;
    m_Vertexes[15].Normal = N4; m_Vertexes[16].Normal = N4; m_Vertexes[17].Normal = N4;
    m_Vertexes[18].Normal = N5; m_Vertexes[19].Normal = N5; m_Vertexes[20].Normal = N5;
    m_Vertexes[21].Normal = N6; m_Vertexes[22].Normal = N6; m_Vertexes[23].Normal = N6;
    m_Vertexes[24].Normal = N7; m_Vertexes[25].Normal = N7; m_Vertexes[26].Normal = N7;
    m_Vertexes[27].Normal = N7; m_Vertexes[28].Normal = N7; m_Vertexes[29].Normal = N7;
    m_Vertexes[30].Normal = N8; m_Vertexes[31].Normal = N8; m_Vertexes[32].Normal = N8;
    m_Vertexes[33].Normal = N8; m_Vertexes[34].Normal = N8; m_Vertexes[35].Normal = N8;

    m_Indexes.resize(36);
    for (unsigned i = 0; i < m_Indexes.size(); ++i)
        m_Indexes[i] = i;
}
