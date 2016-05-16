#include "../stdafx.h"
#include "MeshMathSurface.h"

CMeshMathSurface::CMeshMathSurface()
{
}

CMeshMathSurface::~CMeshMathSurface()
{
}

void CMeshMathSurface::Tesselate()
{
    unsigned long k = 0;
    m_Indexes.resize((m_nVx - 1) * (m_nVy - 1) * 3 * 2);
    for (unsigned long j = 0; j < m_nVy - 1; ++j)
        for (unsigned long i = 0; i < m_nVx - 1; ++i)
        {
            m_Indexes[0 + k] = j * m_nVx + i;
            m_Indexes[1 + k] = (j + 1) * m_nVx + i;
            m_Indexes[2 + k] = j * m_nVx + i + 1;
            m_Indexes[3 + k] = m_Indexes[2 + k];
            m_Indexes[4 + k] = m_Indexes[1 + k];
            m_Indexes[5 + k] = (j + 1) * m_nVx + i + 1;
            k += 6;
        }
}

void CMeshMathSurface::BuildAnalyticSurface(unsigned long nVx, unsigned long nVy, float x0, float y0, float dx, float dy, float(*pFn)(float x, float y), VECTOR4D(*pFnDivergent)(float x, float y, float z))
{
    m_Indexes.clear();
    m_Vertexes.clear();
    m_nVx = nVx;
    m_nVy = nVy;
    float x, y = y0;
    m_Vertexes.resize(m_nVx * m_nVy);
    for (unsigned long j = 0; j < m_nVy; ++j)
    {
        x = x0;
        for (unsigned long i = 0; i < m_nVx; i++)
        {
            // La altura Z se calcula con la funcion que se paso como parametro
            VECTOR4D P{ x, y, pFn(x, y), 1 };
            m_Vertexes[j * m_nVx + i].Position = P;
            m_Vertexes[j * m_nVx + i].Normal = pFnDivergent(x, y, P.z);
            x += dx;
        }
        y += dy;
    }
    Tesselate();
}

void CMeshMathSurface::BuildParametricSurface(unsigned long nVx, unsigned long nVy, float u0, float v0, float du, float dv, VECTOR4D(*pFn)(float u, float v))
{
    m_Indexes.clear();
    m_Vertexes.clear();
    m_nVx = nVx;
    m_nVy = nVy;
    m_Vertexes.resize(m_nVx * m_nVy);
    float u, v = u0;
    for (unsigned long j = 0; j < m_nVy; ++j)
    {
        u = v0;
        for (unsigned long i = 0; i < m_nVx; ++i)
        {
            VECTOR4D P = pFn(u, v);
            m_Vertexes[j * m_nVx + i].Position = P;
            u += du;
        }
        v += dv;
    }
    Tesselate();
}

// Interpolacion bilineal de color
// f(x,y) = Ax + Bxy + Cy + D
void CMeshMathSurface::SetColor(const VECTOR4D &A, const VECTOR4D &B, const VECTOR4D &C, const VECTOR4D &D)
{
    float dcx = 1.0f / m_nVx;
    float dcy = 1.0f / m_nVy;
    float cx, cy = 0.0f;
    for (unsigned long j = 0; j < m_nVy; ++j)
    {
        cx = 0.0f;
        for (unsigned long i = 0; i < m_nVx; ++i)
        {
            m_Vertexes[j * m_nVx + i].Color = Lerp(Lerp(A, B, cx), Lerp(C, D, cx), cy);
            cx += dcx;
        }
        cy += dcy;
    }
}

void CMeshMathSurface::BuildTextureCoords(float u0, float v0, float du, float dv)
{
    float u, v;
    v = v0;
    for (unsigned long j = 0; j < m_nVy; ++j)
    {
        u = u0;
        for (unsigned long i = 0; i < m_nVx; ++i)
        {
            VECTOR4D TexCoord = { u, v, 0, 0 };
            m_Vertexes[j * m_nVx + i].TexCoord = TexCoord;
            u += du;
        }
        v += dv;
    }
}
