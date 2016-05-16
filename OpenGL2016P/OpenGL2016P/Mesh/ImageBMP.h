#pragma once
#include "GL\glew.h"

class CImageBMP
{
    unsigned long m_ulSizeX;    // Pixeles horizontales
    unsigned long m_ulSizeY;    // Pixeles verticales
    struct PIXEL
    {
        unsigned char r, g, b, a;
    };
    PIXEL *m_pBuffer;
protected:
    CImageBMP();
    ~CImageBMP();
public:
    static CImageBMP *CreateBitmapFromFile(char *pszFileName, PIXEL (*pFnAlpha)(PIXEL &P));
    GLuint CreateTexure();
    static void DestroyBitmap(CImageBMP *pBmp);
};
