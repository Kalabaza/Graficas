#include "../stdafx.h"
#include "ImageBMP.h"


CImageBMP::CImageBMP()
{
    m_pBuffer = nullptr;
    m_ulSizeX = 0;
    m_ulSizeY = 0;
}


CImageBMP::~CImageBMP()
{
}


#include <fstream>
#include <Windows.h>
using namespace std;
CImageBMP *CImageBMP::CreateBitmapFromFile(char *pszFileName, PIXEL(*pFnAlpha)(PIXEL &P))
{
    CImageBMP *pNewImage = nullptr;
    fstream file;
    file.open(pszFileName, ios::in | ios::binary);
    if (!file.is_open())
        return nullptr;
    BITMAPFILEHEADER bFh;
    BITMAPINFOHEADER bIh;

    memset(&bFh, 0, sizeof BITMAPFILEHEADER);
    memset(&bIh, 0, sizeof BITMAPINFOHEADER);

    file.read((char*)&bFh.bfType, sizeof(WORD));
    if (bFh.bfType != 'MB')
        return nullptr;
    file.read((char*)&bFh.bfSize, sizeof BITMAPFILEHEADER - sizeof WORD);

    file.read((char*)&bIh.biSize, sizeof DWORD);
    if (bIh.biSize != sizeof BITMAPINFOHEADER)
        return nullptr;
    file.read((char*)&bIh.biWidth, sizeof BITMAPINFOHEADER - sizeof DWORD);

    // Regresa el numero de 32 bytes que se necesitan para leer la imagen
    unsigned long rowLength = 4 * ((bIh.biBitCount * bIh.biWidth + 31) / 32);
    // Ya se tiene la informacion del bitmap en bIh.
    // Se procede a la carga de datos.
    switch(bIh.biBitCount)
    {
    case 1:  // 1bpp - Dicromatico
    {
//        RGBQUAD palette[2];
    }
        break;
    case 4:  // 4bpp - Hexadecacromatico
    {
//        RGBQUAD palette[16];
    }
        break;
    case 8:  // 8bpp - 256 colores
    {
        RGBQUAD palette[256];
        unsigned long colors = bIh.biClrUsed == 0 ? 256 : bIh.biClrUsed;
        file.read((char*)palette, sizeof RGBQUAD * colors);
        unsigned char *pRow = new unsigned char[rowLength];
        pNewImage = new CImageBMP();
        pNewImage->m_ulSizeX = bIh.biWidth;
        pNewImage->m_ulSizeY = bIh.biHeight;
        pNewImage->m_pBuffer = new PIXEL[bIh.biWidth * bIh.biHeight];
        for(long j = 0; j < bIh.biHeight; ++j)
        {
            file.read((char*)pRow, rowLength);
            for(long i = 0; i < bIh.biWidth; ++i)
            {
                PIXEL *p = &pNewImage->m_pBuffer[(bIh.biHeight - j - 1) * bIh.biWidth + i];
                p->b = palette[pRow[i]].rgbBlue;
                p->g = palette[pRow[i]].rgbGreen;
                p->r = palette[pRow[i]].rgbRed;
                p->a = palette[pRow[i]].rgbReserved;
            }
        }
    }
        break;
    case 24: // 24bpp 16M colores
        break;
    case 32: // 32bpp 16M colores + 256 alphas
        break;
    }

    return pNewImage;
}

GLuint CImageBMP::CreateTexure()
{
    GLuint TextID = 0;
    glGenTextures(1, &TextID);
    glBindTexture(GL_TEXTURE_2D, TextID);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, m_ulSizeX, m_ulSizeY, 0, GL_RGBA, GL_UNSIGNED_BYTE, m_pBuffer);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    return TextID;
}

void CImageBMP::DestroyBitmap(CImageBMP *pBmp)
{
    delete[] pBmp->m_pBuffer;
    delete pBmp;
}
