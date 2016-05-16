#pragma once

const int MaxSize = 4;

struct VECTOR4D
{
    union
    {
        struct
        {
            float x, y, z, w;
        };
        struct
        {
            float r, g, b, a;
        };
        float v[MaxSize];
    };
};

struct MATRIX4D
{
    union
    {
        struct
        {
            float m00, m01, m02, m03;
            float m10, m11, m12, m13;
            float m20, m21, m22, m23;
            float m30, m31, m32, m33;
        };
        float m[MaxSize][MaxSize];
        VECTOR4D vec[MaxSize];
        float v[MaxSize*MaxSize];
    };
};

MATRIX4D operator*(const MATRIX4D &A, const MATRIX4D &B);
VECTOR4D operator*(const MATRIX4D &A, const VECTOR4D &V);
VECTOR4D operator*(const VECTOR4D &V, const MATRIX4D &A);
VECTOR4D operator*(const VECTOR4D &A, const VECTOR4D &B);
VECTOR4D operator-(const VECTOR4D &A, const VECTOR4D &B);
VECTOR4D operator+(const VECTOR4D &A, const VECTOR4D &B);
VECTOR4D operator*(const VECTOR4D &A, const float B);
VECTOR4D Cross3(const VECTOR4D &A, const VECTOR4D &B);
float    Dot(const VECTOR4D &A, const VECTOR4D &B);
float    Magnity(const VECTOR4D &A);
VECTOR4D Normalize(const VECTOR4D &A);
MATRIX4D Zero();
MATRIX4D Identity();
MATRIX4D RotationX(float theta);
MATRIX4D RotationY(float theta);
MATRIX4D RotationZ(float theta);
MATRIX4D RotationAxis(float theta, VECTOR4D &Axis);
MATRIX4D Translation(float dx, float dy, float dz);
MATRIX4D Scaling(float sx, float sy, float sz);
MATRIX4D View(VECTOR4D &EyePos, VECTOR4D &Target, VECTOR4D &Up);
MATRIX4D PerspectiveWidthHeightRH(float nWidth, float nHeight, float zNear, float zFar);
MATRIX4D Transpose(const MATRIX4D &M);
MATRIX4D FastInverse(const MATRIX4D &M);
VECTOR4D Lerp(const VECTOR4D &A, const VECTOR4D &B, float u);   // I = A + u * (B - A)
bool PtInTriangle(const VECTOR4D &V0, const VECTOR4D &V1, const VECTOR4D &V2, const VECTOR4D &P);
bool PtInTriangleBarycentric(const VECTOR4D &V0, const VECTOR4D &V1, const VECTOR4D &V2, const VECTOR4D &P, float &w0, float &w1, float &w2);
float Inverse(MATRIX4D& M, MATRIX4D& R);
void PlaneIntersect(VECTOR4D &RayOrigin, VECTOR4D &RayDir, VECTOR4D &Plane, float &n, float &d); // Source, Direction, Plane eq. u = n/d, V'=V0+u*RayDir
bool RayCastOnTriangle(VECTOR4D &V0, VECTOR4D &V1, VECTOR4D &V2, VECTOR4D &RayOrigin, VECTOR4D &RayDir, VECTOR4D &Intersection);
void BuildRayFromPerspective(MATRIX4D &PV, float x, float y, VECTOR4D &RayOrigin, VECTOR4D &RayDir);
MATRIX4D Orthogonalize(MATRIX4D &M);
