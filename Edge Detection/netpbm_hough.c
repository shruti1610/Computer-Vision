/*
 * netpbm_hough.c
 *
 *  Created on: Nov 11, 2019
 *      Author: shrutighelani
 */


// netpbm_hough.c
// Test and demo program for the Hough transform using the netpbm.c library.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "netpbm.h"

#define MIN(X,Y) ((X)<(Y)?(X):(Y))
#define MAX(X,Y) ((X)>(Y)?(X):(Y))
#define PI 3.141592653589793


// Linearly convolve the intensity information of an image with a convolution filter (both
// stored in Matrices). The center of the filter serves as its anchor for storing the result.
// For even number of rows or columns, the anchor will be rounded to the left or upward, resp.
// The output is a matrix of the same size as the original image, with zeros filling the
// border regions that are unreachable by the filter.
Matrix convolve(Matrix img, Matrix filter)
{
	int m, n, k, l, mOffset = filter.height / 2, nOffset = filter.width / 2;
	double sum;
	Matrix result = createMatrix(img.height, img.width);

	for (m = filter.height - 1; m < img.height; m++)
		for (n = filter.width - 1; n < img.width; n++)
		{
			sum = 0.0;
			for (k = 0; k < filter.width; k++)
				for (l = 0; l < filter.height; l++)
					sum += img.map[m - k][n - l] * filter.map[k][l];
			result.map[m - mOffset][n - nOffset] = sum;
		}
	return result;
}


// Build a Hough parameter map for matrix mxSpatial for detecting straight lines.
// Rows in this map represent the normal alpha and columns represent the distance d from the origin.
// Increasing the size of the map in each dimension improves the resolution of the corresponding parameter.


// FOR CIRCLE

Matrix houghTransformLines(Matrix mxSpatial, int mapHeight, int mapWidth, Matrix output, int hori, int vert)
{
	int m, n, angle, verticalDist, horizontalDist;
	double alpha;
	Matrix mxParam = createMatrix(mxSpatial.height,mxSpatial.width);
	Matrix sincos = createMatrix(mapHeight, 2);

	// Generate lookup table for sin and cos values to speed up subsequent computation.
	for (angle = 0; angle < mapHeight; angle++)
	{
		/*alpha = -0.5*PI + 1.5*PI*(double) angle/(double) mapHeight;
		sincos.map[angle][0] = (double) mapWidth/maxD*sin(alpha);
		sincos.map[angle][1] = (double) mapWidth/maxD*cos(alpha);*/

		alpha = angle * PI/180;
		sincos.map[angle][0] = (double) sin(alpha);
		sincos.map[angle][1] = (double) cos(alpha);
	}

	for (m = 0; m < mxSpatial.height; m++)
		for (n = 0; n < mxSpatial.width; n++)
			/*for (angle = 0; angle < mapHeight; angle++)
			{
				dist = (int) ((double) m*sincos.map[angle][0] + (double) n*sincos.map[angle][1] + 0.5);
				if (dist >= 1)
					mxParam.map[angle][dist] += mxSpatial.map[m][n];
			}

	deleteMatrix(sincos);
	return mxParam;
}*/
			if(mxSpatial.map[m][n] >= 200)
			{

					for (angle =0; angle < mapHeight; angle++)
					{
						verticalDist =(int) (m - (double)vert*sincos.map[angle][0]);
						horizontalDist = (int)(n - (double)hori*sincos.map[angle][1]);
						if(verticalDist>0 && verticalDist<mxSpatial.height && horizontalDist>0 &&  horizontalDist<mxSpatial.width)
						{
							mxParam.map[verticalDist][horizontalDist] += mxSpatial.map[m][n];
							output.map[verticalDist][horizontalDist] +=mxSpatial.map[m][n];
						}
					}
			}

			deleteMatrix(sincos);
			return mxParam;
	}


// Test whether entry (m, n) in matrix mx is a local maximum, i.e., is not exceeded by any of its
// maximally 8 neighbors. Return 1 if true, 0 otherwise.
int isLocalMaximum(Matrix mx, int m, int n)
{
	double strength = mx.map[m][n];
	int i, j;
	int iMin = (m == 0)? 0:(m - 1);
	int iMax = (m == mx.height -1)? m:(m + 1);
	int jMin = (n == 0)? 0:(n - 1);
	int jMax = (n == mx.width -1)? n:(n + 1);

	for (i = iMin; i <= iMax; i++)
		for (j = jMin; j <= jMax; j++)
			if (mx.map[i][j] > strength)
				return 0;
	return 1;
}


// Insert a new entry, consisting of vPos, hPos, and strength, into the list of maxima mx.
void insertMaxEntry(Matrix mx, int vPos, int hPos, double strength)
{
	int m, n = mx.width - 1;

	while (n > 0 && mx.map[2][n - 1] < strength)
	{
		for (m = 0; m < 3; m++)
			mx.map[m][n] = mx.map[m][n - 1];
		n--;
	}
	mx.map[0][n] = (double) vPos;
	mx.map[1][n] = (double) hPos;
	mx.map[2][n] = strength;
}


// Delete entry number i from the list of maxima mx.
void deleteMaxEntry(Matrix mx, int i)
{
	int m, n;

	for (n = i; n < mx.width - 1; n++)
		for (m = 0; m < 3; m++)
			mx.map[m][n] = mx.map[m][n + 1];

	mx.map[2][mx.width - 1] = -1.0;
}


// Find the <number> highest maxima in a Hough parameter map that are separated by a Euclidean distance
// of at least <minSeparation> in the map. The result is a three-row matrix with each column representing
// the row, the column, and the strength of one maximum, in descending order of strength.
Matrix findHoughMaxima(Matrix mx, int number, double minSeparation)
{
	int j, m, n, k, l, r, index, do_not_insert;
	double minSepSquare = SQR(minSeparation);
	double strength;
	Matrix maxima = createMatrix(3, number);

	for (j = 0; j < number; j++)
		maxima.map[2][j] = -1.0;

	for (m = 0; m < mx.height; m++)
		for (n = 0; n < mx.width; n++)
		{
			strength = mx.map[m][n];
			index = 0;
			do_not_insert = 0;
			if (strength > 0.0 && strength > maxima.map[2][number - 1] && isLocalMaximum(mx, m, n))
			{
				while (index < number && !do_not_insert)
				{
					if (maxima.map[2][index] > 0.0 && SQR(maxima.map[0][index] - (double) m) + SQR(maxima.map[1][index] - (double) n) < minSepSquare)
						if (strength > maxima.map[2][index])
							deleteMaxEntry(maxima, index);
						else
							do_not_insert = 1;
					index++;
				}
				if (!do_not_insert)
					insertMaxEntry(maxima, m, n, strength);
				index = 0;
			}
		}
	return maxima;
}


// Read image "desk.ppm" and write Hough transform related output images.
void main()
{
	int i, m, n, m1, n1, m2, n2,verticalDist, horizontalDist;
	double hFilter[3][3] = { { 1.0, 0.0, -1.0 },{ 2.0, 0.0, -2.0 },{ 1.0, 0.0, -1.0 } };
	double vFilter[3][3] = { { 1.0, 2.0, 1.0 },{ 0.0, 0.0, 0.0 },{ -1.0, -2.0, -1.0 } };
	Matrix hSobel = createMatrixFromArray(&hFilter[0][0], 3, 3);
	Matrix vSobel = createMatrixFromArray(&vFilter[0][0], 3, 3);
	Image inputImage = readImage("desk.ppm");
	Matrix inputMatrix = image2Matrix(inputImage);
	Matrix hEdgeMatrix = convolve(inputMatrix, hSobel);
	Matrix vEdgeMatrix = convolve(inputMatrix, vSobel);
	Image edgeImage, houghImage;
	Matrix edgeMatrix = createMatrix(inputImage.height, inputImage.width);
	Matrix houghMatrix, maxMatrix;
	double maxLength, alpha, dist;
	int hr,vr;

	houghMatrix = createMatrix(inputImage.height, inputImage.width);

	int numLines = 5, alphaBins = 360, dBins = 500;
	double minPeakDistance = 50.0;

	for (m = 0; m < edgeMatrix.height; m++)
		for (n = 0; n < edgeMatrix.width; n++)
			edgeMatrix.map[m][n] = sqrt(SQR(hEdgeMatrix.map[m][n]) + SQR(vEdgeMatrix.map[m][n]));

	edgeImage = matrix2Image(edgeMatrix, 1, 1.0);
	writeImage(edgeImage, "Desk_edge.pgm");

	hr =edgeImage.width/2;
	vr=	edgeImage.height/2;

	Matrix houghArray;

	houghArray = houghTransformLines(edgeMatrix, alphaBins, dBins, houghMatrix,hr,vr);
	houghImage = matrix2Image(houghMatrix, 1, 1.0);
	writeImage(houghImage, "Desk_hough.pgm");
	printf("initiated maxima\n");

	maxMatrix = findHoughMaxima(houghArray, numLines, minPeakDistance);

	printf("ellipse started\n");
	for (i = 0; i < numLines; i++)
		ellipse(houghImage, maxMatrix.map[0][i], maxMatrix.map[1][i],20,20, 2, 10, 7, 255, 255, 255, 0);
	writeImage(houghImage, "Desk_hough_max.ppm");

	printf("completed ellipse\n");



	/*for (i = 0; i < numLines; i++)
	{
		alpha = -0.5*PI + 1.5*PI*maxMatrix.map[0][i]/(double) houghMatrix.height;
		dist = maxLength*maxMatrix.map[1][i]/(double) houghMatrix.width;
		m1 = (int) (dist*sin(alpha) - maxLength*cos(alpha) + 0.5);
		n1 = (int) (dist*cos(alpha) + maxLength*sin(alpha) + 0.5);
		m2 = (int) (dist*sin(alpha) + maxLength*cos(alpha) + 0.5);
		n2 = (int) (dist*cos(alpha) - maxLength*sin(alpha) + 0.5);
		line(inputImage, m1, n1, m2, n2, 2, 18, 10, 30, 30, 30, 0);
	}	writeImage(inputImage, "desk_hough_lines.ppm");*/

	for (i = 0; i < numLines; i++)
	{
		/*alpha = maxMatrix.map[0][i]* PI/180 ;
		verticalDist = maxLength*maxMatrix.map[1][i]/(double) houghMatrix.width;
		horizontalDist=maxLength*maxMatrix.map[1][i]/(double) houghMatrix.width;
		//m1 = (int) (dist*sin(alpha) - maxLength*cos(alpha));
		//n1 = (int) (dist*cos(alpha) + maxLength*sin(alpha));
		//m2 = (int) (dist*sin(alpha) + maxLength*cos(alpha));
		//n2 = (int) (dist*cos(alpha) - maxLength*sin(alpha));
		m1 =(int) (double)verticalDist*sin(alpha); // sin vertical
		n1 = (int)(double)horizontalDist*cos(alpha);*/
		ellipse(inputImage,maxMatrix.map[0][i] ,maxMatrix.map[1][i] ,20, 20,2,10,7,255,255,255,0);





	}writeImage(inputImage, "Desk_hough_lines.ppm");


	deleteMatrix(edgeMatrix);
	deleteMatrix(houghMatrix);
	deleteImage(inputImage);
	deleteImage(edgeImage);
	deleteImage(houghImage);
}
