#pragma once

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
        float v[4];
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
        float m[4][4];
        VECTOR4D vec[4];
        float v[16];
    };
};

MATRIX4D operator*(const MATRIX4D& A, const MATRIX4D &B);
VECTOR4D operator*(const MATRIX4D& A, const VECTOR4D& V);
VECTOR4D operator*(const VECTOR4D& V, const MATRIX4D& A);
VECTOR4D operator*(const VECTOR4D& A, const VECTOR4D& B);
VECTOR4D operator-(const VECTOR4D& A, const VECTOR4D &B);
VECTOR4D operator+(const VECTOR4D& A, const VECTOR4D &B);
VECTOR4D Cross3(const VECTOR4D&A, const VECTOR4D &B);
float    Dot(const VECTOR4D& A, const VECTOR4D& B);
float    Magnity(const VECTOR4D& A);
VECTOR4D Normalize(const VECTOR4D& A);
MATRIX4D Zero();
MATRIX4D Identity();
MATRIX4D RotationX(float theta);
MATRIX4D RotationY(float theta);
MATRIX4D RotationZ(float theta);
MATRIX4D Translation(float dx, float dy, float dz);
MATRIX4D Scaling(float sx, float sy, float sz);
MATRIX4D View(VECTOR4D& EyePos,
VECTOR4D& Target, VECTOR4D& Up);
