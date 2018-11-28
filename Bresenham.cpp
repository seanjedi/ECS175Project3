#include "Bresenham.h"
#include <math.h>
#include <iostream>

using namespace std;
void makePixel(int x, int y, float* PixelBuffer, float intensity);
inline void max(int& a, int& b, float& c, float& d) {
	if (a > b) {
		int temp = a;
		a = b;
		b = temp;
		float temp2 = c;
		c = d;
		d = temp2;
	}
}

void Bresenham(int x1, int x2, int y1, int y2, float* PixelBuffer, int windowSizeX, float light1, float light2) {
	float intensity;
	if (x1 == x2) {
		max(y1, y2, light1, light2);
		for (int i = y1; i < y2; i++) {;
			intensity = ((float(y2) - float(i)) / (float(y2) - float(y1)) * light1) 
				+ ((float(i) - float(y1)) / (float(y2) - float(y1)) * light2);
			makePixel(x1, i, PixelBuffer, intensity);
		}
	}
	else if (y1 == y2) {
		max(x1, x2, light1, light2);
		for (int i = x1; i < x2; i++) {
			intensity = ((float(x2) - float(i)) / (float(x2) - float(x1)) * light1) 
				+ ((float(i) - float(x1)) / (float(x2) - float(x1)) * light2);
			makePixel(i, y1,  PixelBuffer, intensity);
		}
	}
	else {
		int dx = fabs(x2 - x1), dy = fabs(y2 - y1);
		int x, y;

		float m = (float(y2 - y1) / float(x2 - x1));
		if (m >= 1) {
			if (y1 > y2) {
				x = x2;
				y = y2;
				max(y1, y2, light1, light2);
			}
			else {
				x = x1;
				y = y1;
			}
			int Dx2 = 2 * dx, Dx2y = 2 * (dx - dy);
			int p = 2 * dx - dy;
			intensity = ((float(y2) - float(y)) / (float(y2) - float(y1)) * light1) 
				+ ((float(y) - float(y1)) / (float(y2) - float(y1)) * light2);
			makePixel(x, y, PixelBuffer, intensity);
			while (y < y2) {
				y++;
				if (p < 0)
					p += Dx2;
				else {
					x++;
					p += Dx2y;
				}
				intensity = ((float(y2) - float(y)) / (float(y2) - float(y1)) * light1) 
					+ ((float(y) - float(y1)) / (float(y2) - float(y1)) * light2);
				makePixel(x, y, PixelBuffer, intensity);
			}
		}
		else if (m > 0 && m < 1) {
			if (x1 > x2) {
				x = x2;
				y = y2;
				max(x1, x2, light1, light2);
			}
			else {
				x = x1;
				y = y1;
			}
			int Dy2 = 2 * dy, Dy2x = 2 * (dy - dx);
			int p = 2 * dy - dx;
			intensity = ((float(x2) - float(x)) / (float(x2) - float(x1)) * light1) 
				+ ((float(x) - float(x1)) / (float(x2) - float(x1)) * light2);
			makePixel(x, y, PixelBuffer, intensity);
			while (x < x2) {
				x++;
				if (p < 0)
					p += Dy2;
				else {
					y++;
					p += Dy2x;
				}
				intensity = ((float(x2) - float(x)) / (float(x2) - float(x1)) * light1) 
					+ ((float(x) - float(x1)) / (float(x2) - float(x1)) * light2);
				makePixel(x, y, PixelBuffer, intensity);
			}
		}
		else if (m <= -1) {
			if (y1 > y2) {
				x = x2;
				y = y2;
				max(y1, y2, light1, light2);
			}
			else {
				x = x1;
				y = y1;
			}
			int Dx2 = 2 * dx, Dx2y = 2 * (dx - dy);
			int p = 2 * dx - dy;
			intensity = ((float(y2) - float(y)) / (float(y2) - float(y1)) * light1) 
				+ ((float(y) - float(y1)) / (float(y2) - float(y1)) * light2);
			makePixel(x, y, PixelBuffer, intensity);
			while (y < y2) {
				y++;
				if (p < 0)
					p += Dx2;
				else {
					x--;
					p += Dx2y;
				}
				intensity = ((float(y2) - float(y)) / (float(y2) - float(y1)) * light1) 
					+ ((float(y) - float(y1)) / (float(y2) - float(y1)) * light2);
				makePixel(x, y, PixelBuffer, intensity);
			}
		}
		else if (m < 0 && m > -1) {
			if (x1 > x2) {
				x = x2;
				y = y2;
				max(x1, x2, light1, light2);
			}
			else {
				x = x1;
				y = y1;
			}
			int Dy2 = 2 * dy, Dy2x = 2 * (dy - dx);
			int p = 2 * dy - dx;
			intensity = ((float(x2) - float(x)) / (float(x2) - float(x1)) * light1) 
				+ ((float(x) - float(x1)) / (float(x2) - float(x1)) * light2);
			makePixel(x, y, PixelBuffer, intensity);
			while (x < x2) {
				x++;
				if (p < 0)
					p += Dy2;
				else {
					y--;
					p += Dy2x;
				}
				intensity = ((float(x2) - float(x)) / (float(x2) - float(x1)) * light1) 
					+ ((float(x) - float(x1)) / (float(x2) - float(x1)) * light2);
				makePixel(x, y, PixelBuffer, intensity);
			}
		}
	}
}