#pragma once

#include "Mesh.h"

class CMeshMathSurface : public CMesh
{
protected:
    unsigned long m_nVx;
    unsigned long m_nVy;

    void Tesselate();               // Funcion que genera todos los indices necesarios para generar una figura geometrica
public:
    CMeshMathSurface();
    ~CMeshMathSurface();
    // Funcion para generar superficies por medio de funciones geometricas
    void BuildAnalyticSurface(unsigned long nVx, unsigned long nVy, // Resolucion de la malla
                              float x0, float y0,                   // x0 X y0 = z0
                              float dx, float dy,                   // Diferencial de la distancia entre los vertices
                              float(*pFn)(float x, float y),        // Apuntador a la funcion que genera los valores de Z
                              VECTOR4D(*pFnDivergent)(float x, float y, float z));
    // Funcion que puede producir superficies cerradas
    void BuildParametricSurface(unsigned long nVx, unsigned long nVy, // Resolucion de la malla
                                float u0, float v0,
                                float du, float dv,                   // Diferencial de la distancia entre los vertices
                                VECTOR4D(*pFn)(float u, float v));
    void BuildTextureCoords(float u0, float v0, float du, float dv);
    void SetColor(const VECTOR4D &A, const VECTOR4D &B, const VECTOR4D &C, const VECTOR4D &D);
};
