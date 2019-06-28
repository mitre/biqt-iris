// #######################################################################
// NOTICE
//
// This software (or technical data) was produced for the U.S. Government
// under contract, and is subject to the Rights in Data-General Clause
// 52.227-14, Alt. IV (DEC 2007).
//
// Copyright 2019 The MITRE Corporation. All Rights Reserved.
// #######################################################################
//
//        FILE: ImageOps.cpp
// DESCRIPTION: Class that provides open source Defocus
//              and Overall Contrast quality metrics for
//              .bmp iris images.
//        DATE: 13 June 2012
//
//---------------------------------------------------------------------

#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <math.h>
#include <stdexcept>
#include <limits>

#include "ImageOps.h"
#ifdef _WIN32
#include <windows.h>
#endif

/* static */ const int
    MFilter::FocusKernel[FOCUS_KERNEL_DIM][FOCUS_KERNEL_DIM] = {
        {-1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1},
        {-1, -1, 3, 3, 3, 3, -1, -1},     {-1, -1, 3, 3, 3, 3, -1, -1},
        {-1, -1, 3, 3, 3, 3, -1, -1},     {-1, -1, 3, 3, 3, 3, -1, -1},
        {-1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1}};
/* static */ const int
    MFilter::HalfFocusKernel[HALF_FOCUS_KERNEL_DIM][HALF_FOCUS_KERNEL_DIM] = {
        {-1, -1, -1, -1}, {-1, 3, 3, -1}, {-1, 3, 3, -1}, {-1, -1, -1, -1}};
/* static */ int MFilter::PolarSegMax = MFilter::maxsegs;
/* static */ int MFilter::PLowerSegTrial = 63;
/* static */ int MFilter::PUpperSegTrial = 35;
/* static */ double MFilter::dMaxR = (double)MFilter::MaxR;
/* static */ int MFilter::seg90 = 0;
/* static */ double MFilter::singlesegment = 0;
/* static */ int MFilter::LowerSegTrial = 47;
/* static */ int MFilter::UpperSegTrial = 35;

/**
 *
 */
MFilter::MFilter()
    : MathTablesDefined(false), DefocusScore(0), NormDefocusScore(0),
      DefocusScoreRef(0), IrisCenterX(0), IrisCenterY(0), IrisRadius(0),
      RoughIrisRadius(0), RoughIriscx(0), RoughIriscy(0), RoughPcx(320),
      RoughPcy(240), RoughPrad(50), RoughPdiam(100), bestresponse(0),
      rightresponse(0), leftresponse(0), maxIrisFindLeftEdgeResponse(0),
      bestdarkiriscenter(0), maxIrisFindRightEdgeResponse(0),
      MaxPupilEdgeResponse(0), edgeresponse(0), maxVertEdgeResponse(0),
      leftedgeIrisCenterx(0), leftedgeIrisCentery(0), leftedgeIrisDiameter(0),
      rightedgeIrisCenterx(0), rightedgeIrisCentery(0),
      rightedgeIrisDiameter(0), UsableIrisAreaPercent(0),
      ContrastLowLimit(30.0), ContrastUpperLimit(50.0),
      NContrastValueAtLowLimit(0), NContrastValueAtUpperLimit(1),
      NContrastFactor(1), DefocusLowLimit(60), DefocusUpperLimit(80),
      NDefocusValueAtLowLimit(0), NDefocusValueAtUpperLimit(1),
      NDefocusFactor(1), IDLowLimit(170), IDMedLimit1(200), IDMedLimit2(350),
      IDUpperLimit(400), NIDiamFactor(1), ISGSMeanLowLimit(20),
      ISGSMeanUpperLimit(30), NISGSValueAtLowLimit(0),
      NISGSValueAtUpperLimit(1), NISGSFactor(1), IrisVisLowLimit(4),
      IrisVisUpperLimit(7), IPGSDiffLowLimit(15), IPGSDiffUpperLimit(30),
      NIPGSDiffFactor(1), NPupilIrisRatio(1),
      PupilIrisRatioLowLimit((float)0.20), PupilIrisRatioHighLimit((float)0.60),
      NPupilIrisRatioFactor(1), NIrisVisFactor(1), CombinedQualityLowLimit(0),
      CombinedQualityUpperLimit((float)0.80), OverallQuality(0), MinRtrial(28),
      MaxRtrial(42), MinPrtrial(8), MaxPrtrial(32), pupilCx(0), pupilCy(0),
      pupilRad(5), pupilD(10), ContrastScore(0), AverageImageIntensity(0),
      LMargin(0), RMargin(0), TMargin(0), BMargin(0),
      LRMarginMetricTotalDeduct(0.0), TBMarginMetricTotalDeduct(0.0),
      OverallMargin(1.0), ISGSDiffMeanAvg(0.0), ISGSDiffMeanMax(0.0),
      IrisPupilGSDiff(0), pointresponse(0), maxpointresponse(2000)
{
    uppersegpointresponses = new int[LowerSegTrial];
    maxuppersegresponses = new int[LowerSegTrial];
    Framebin = new unsigned char[maxwidth * maxheight];
    ISGSVals = new int[256];
    PupilVals = new int[256];
    irispts = new Point[101];
    cos = new int[101];

    Rawbin = new unsigned char *[maxwidth];
    ValF16Edge = new int *[maxwidth];
    ValF16EdgeBytes = new unsigned char *[maxwidth];
    TransposedRawbin = new unsigned char *[maxwidth];
    MaskedRawbin = new unsigned char *[maxwidth];
    Rawvals = new int *[maxwidth];

    int i;

    for (i = 0; i < maxwidth; i++) {
        Rawbin[i] = new unsigned char[maxheight];
        ValF16Edge[i] = new int[maxheight];
        ValF16EdgeBytes[i] = new unsigned char[maxheight];
        TransposedRawbin[i] = new unsigned char[maxheight];
        MaskedRawbin[i] = new unsigned char[maxheight];
        Rawvals[i] = new int[maxheight];
    }

    for (i = 0; i < maxwidth; i++) {
        memset(Rawbin[i], 0, sizeof(unsigned char) * maxheight);
        memset(ValF16Edge[i], 0, sizeof(int) * maxheight);
        memset(ValF16EdgeBytes[i], 0, sizeof(unsigned char) * maxheight);
        memset(TransposedRawbin[i], 0, sizeof(unsigned char) * maxheight);
        memset(MaskedRawbin[i], 0, sizeof(unsigned char) * maxheight);
        memset(Rawvals[i], 0, sizeof(int) * maxheight);
    }

    DSRawbin = new unsigned char *[(int)(maxwidth / 4)];
    DSvalbin = new int *[(int)(maxwidth / 4)];
    NegDSRawbin = new unsigned char *[(int)(maxwidth / 4)];
    NegDSvalbin = new int *[(int)(maxwidth / 4)];
    Edgevalbin = new int *[(int)(maxwidth / 4)];
    PupilEdgebin = new int *[(int)(maxwidth / 4)];
    TransposedEdgebin = new int *[(int)(maxwidth / 4)];
    EdgeRawbin = new unsigned char *[(int)(maxwidth / 4)];
    VertEdge = new int *[(int)(maxwidth / 4)];
    VertEdgebin = new unsigned char *[(int)(maxwidth / 4)];
    for (i = 0; i < (int)(maxwidth / 4); i++) {
        DSRawbin[i] = new unsigned char[(int)(maxheight / 4)];
        DSvalbin[i] = new int[(int)(maxheight / 4)];
        NegDSRawbin[i] = new unsigned char[(int)(maxheight / 4)];
        NegDSvalbin[i] = new int[(int)(maxheight / 4)];
        Edgevalbin[i] = new int[(int)(maxheight / 4)];
        PupilEdgebin[i] = new int[(int)(maxheight / 4)];
        TransposedEdgebin[i] = new int[(int)(maxheight / 4)];
        EdgeRawbin[i] = new unsigned char[(int)(maxheight / 4)];
        VertEdge[i] = new int[(int)(maxheight / 4)];
        VertEdgebin[i] = new unsigned char[(int)(maxheight / 4)];
    }

    for (i = 0; i < (int)(maxwidth / 4); i++) {
        memset(DSRawbin[i], 0, sizeof(unsigned char) * maxheight / 4);
        memset(DSvalbin[i], 0, sizeof(int) * maxheight / 4);
        memset(NegDSRawbin[i], 0, sizeof(unsigned char) * maxheight / 4);
        memset(Edgevalbin[i], 0, sizeof(int) * maxheight / 4);
        memset(PupilEdgebin[i], 0, sizeof(int) * maxheight / 4);
        memset(TransposedEdgebin[i], 0, sizeof(unsigned char) * maxheight / 4);
        memset(EdgeRawbin[i], 0, sizeof(unsigned char) * maxheight / 4);
        memset(VertEdge[i], 0, sizeof(int) * maxheight / 4);
        memset(VertEdgebin[i], 0, sizeof(unsigned char) * maxheight / 4);
    }

    int j;
    trigvals = new int **[maxsegs + 1];
    radials = new int **[maxsegs + 1];
    for (i = 0; i < maxsegs + 1; i++) {
        trigvals[i] = new int *[MaxR + 1];
        radials[i] = new int *[MaxR + 1];
        for (j = 0; j < MaxR + 1; j++) {
            trigvals[i][j] = new int[2];
            radials[i][j] = new int[2];
        }
    }

    for (i = 0; i < maxsegs + 1; i++) {
        for (j = 0; j < MaxR + 1; j++) {
            memset(trigvals[i][j], 0, sizeof(int) * 2);
            memset(radials[i][j], 0, sizeof(int) * 2);
        }
    }
}

/**
 *
 */
MFilter::~MFilter()
{
    delete[] uppersegpointresponses;
    delete[] maxuppersegresponses;
    delete[] Framebin;
    delete[] ISGSVals;
    delete[] PupilVals;
    delete[] irispts;
    delete[] cos;

    int i;

    for (i = 0; i < maxwidth; i++) {
        delete[] Rawbin[i];
        delete[] ValF16Edge[i];
        delete[] ValF16EdgeBytes[i];
        delete[] TransposedRawbin[i];
        delete[] MaskedRawbin[i];
        delete[] Rawvals[i];
    }
    delete[] Rawbin;
    delete[] ValF16Edge;
    delete[] ValF16EdgeBytes;
    delete[] TransposedRawbin;
    delete[] MaskedRawbin;
    delete[] Rawvals;

    for (i = 0; i < (int)(maxwidth / 4); i++) {
        delete[] DSRawbin[i];
        delete[] DSvalbin[i];
        delete[] NegDSRawbin[i];
        delete[] NegDSvalbin[i];
        delete[] Edgevalbin[i];
        delete[] PupilEdgebin[i];
        delete[] TransposedEdgebin[i];
        delete[] EdgeRawbin[i];
        delete[] VertEdge[i];
        delete[] VertEdgebin[i];
    }
    delete[] DSRawbin;
    delete[] DSvalbin;
    delete[] NegDSRawbin;
    delete[] NegDSvalbin;
    delete[] Edgevalbin;
    delete[] PupilEdgebin;
    delete[] TransposedEdgebin;
    delete[] EdgeRawbin;
    delete[] VertEdge;
    delete[] VertEdgebin;

    int j;
    for (i = 0; i < maxsegs + 1; i++) {
        for (j = 0; j < MaxR + 1; j++) {
            delete[] trigvals[i][j];
            delete[] radials[i][j];
        }
        delete[] trigvals[i];
        delete[] radials[i];
    }
    delete[] trigvals;
    delete[] radials;
}

/**
 *
 * @param rawimg
 * @param imgwidth
 * @param imgheight
 * @return
 */
int MFilter::Defocus(unsigned char **rawimg, int imgwidth, int imgheight)
{
    try {
        int kwidth = (int)FOCUS_KERNEL_DIM;
        int kheight = (int)FOCUS_KERNEL_DIM;
        int marwid = kwidth / 2;
        int marheight = kheight / 2;
        int ptresponse = 0;
        int totalresponse = 0;
        int trialpt = 0;
        int numpixels = imgwidth * imgheight;
        int numpoints = 0;

        if (imgwidth == 0 || imgheight == 0) {
            DefocusScore = 0;
            NormDefocusScore = 0;
            return ERROR_ZERO_IMAGESIZE;
        }

        if (imgwidth > maxwidth || imgheight > maxheight) {
            DefocusScore = 0;
            NormDefocusScore = 0;
            return ERROR_IMAGESIZE_MAX_EXCEEDED;
        }

        int idx1;
        int idy1;
        int jdx;
        int jdy;
        for (idx1 = marwid * 4; idx1 < imgwidth - marheight * 4;
             idx1 += 4) // idx1++, idx1++, idx1++, idx1++)
        {
            for (idy1 = marheight * 4; idy1 < imgheight - marheight * 4;
                 idy1 += 4) // idy1++, idy1++, idy1++, idy1++)
            {
                // operate the kernel at the selected trial point
                ptresponse = 0;
                for (jdx = 0; jdx < kwidth; jdx++) {
                    for (jdy = 0; jdy < kheight; jdy++) {
                        trialpt = (int)
                            rawimg[idx1 - marwid + jdx][idy1 - marheight + jdy];
                        ptresponse =
                            (trialpt * FocusKernel[jdx][jdy]) + ptresponse;
                    }
                }
                totalresponse += (int)abs(ptresponse);
                numpoints++;
            }
        }
        DefocusScore = (int)(totalresponse / numpoints);
        NormDefocusScore =
            (float)((float)(DefocusScore * DefocusScore) /
                    (float)(DefocusScore * DefocusScore + 90000.0));
    }
    catch (std::exception &ex) {
        DefocusScore = 0;
        NormDefocusScore = 0;
        std::cerr << "exception caught: " << ex.what() << std::endl;
        return ERROR_UNSPECIFIED_ERROR;
    }
    return DefocusScore;
}

/**
 *
 * @param framebytes
 * @param width
 * @param height
 * @return
 */
int MFilter::Defocus(unsigned char *framebytes, int width, int height)
{
    unsigned char **MyTwoDArray;
    MyTwoDArray = new unsigned char *[width];
    int i;
    for (i = 0; i < width; i++) {
        MyTwoDArray[i] = new unsigned char[height];
    }
    RawByteFrameToTwoDimArray(framebytes, MyTwoDArray, width, height);
    int result = Defocus(MyTwoDArray, width, height);
    for (i = 0; i < width; i++) {
        delete[] MyTwoDArray[i];
    }
    delete[] MyTwoDArray;
    return result;
}

/**
 *
 * @param rawimg2
 * @param imgwidth2
 * @param imgheight2
 * @return
 */
int MFilter::Contrast(unsigned char **rawimg2, int imgwidth2, int imgheight2)
{
    int ndx4;
    int ndy4;
    double totalCSquare = 0;
    double totalC = 0;
    double avgC = 0;
    double vari;

    try {
        double numpix = (double)imgwidth2 * (double)imgheight2;
        if ((imgwidth2 > maxwidth) || (imgheight2 > maxheight)) {
            return ERROR_IMAGESIZE_MAX_EXCEEDED;
        }
        for (ndy4 = 0; ndy4 < imgheight2; ndy4++) {
            for (ndx4 = 0; ndx4 < imgwidth2; ndx4++) {
                totalC += (double)rawimg2[ndx4][ndy4];
            }
        }
        if (numpix == 0) {
            return ERROR_ZERO_IMAGESIZE;
        }
        avgC = totalC / numpix;

        for (ndy4 = 0; ndy4 < imgheight2; ndy4++) {
            for (ndx4 = 0; ndx4 < imgwidth2; ndx4++) {
                vari = (double)rawimg2[ndx4][ndy4] - avgC;
                totalCSquare += vari * vari;
            }
        }
        ContrastScore = (int)sqrt(totalCSquare / numpix);
        AverageImageIntensity = (int)avgC;
    }
    catch (std::exception &ex) {
        ContrastScore = 0;
        std::cerr << "exception caught: " << ex.what() << std::endl;
        return ERROR_UNSPECIFIED_ERROR;
    }
    return ContrastScore;
}

/**
 *
 * @param framebytes
 * @param imgwidth
 * @param imgheight
 * @return
 */
int MFilter::Contrast(unsigned char *framebytes, int imgwidth, int imgheight)
{
    unsigned char **MyTwoDArray;
    MyTwoDArray = new unsigned char *[imgwidth];
    int i;
    for (i = 0; i < imgwidth; i++) {
        MyTwoDArray[i] = new unsigned char[imgheight];
    }
    RawByteFrameToTwoDimArray(framebytes, MyTwoDArray, imgwidth, imgheight);
    int result = Contrast(MyTwoDArray, imgwidth, imgheight);
    for (i = 0; i < imgwidth; i++) {
        delete[] MyTwoDArray[i];
    }
    delete[] MyTwoDArray;
    return result;
}

/**
 *
 * @param singlearray
 * @param twodimarray
 * @param width
 * @param height
 */
void MFilter::RawByteFrameToTwoDimArray(const unsigned char *singlearray,
                                        unsigned char **twodimarray, int width,
                                        int height)
{
    try {
        int xdx;
        int modx;
        int mody;
        int cntr = 0;

        for (int ydx = 0; ydx < height; ydx++) {
            for (xdx = 0; xdx < width; xdx++) {
                mody = cntr / width;
                modx = cntr % width; // modulo division to get the
                                     // cntr-(trunc(cntr/width))*width
                twodimarray[xdx][ydx] = singlearray[cntr];
                cntr++;
            }
        }
    }
    catch (std::exception &ex) {
        std::cerr << "exception caught: " << ex.what() << std::endl;
    }
}

/**
 *
 * @param imgwidth1
 * @param imgheight1
 * @param IcX
 * @param IcY
 * @param ID
 */
void MFilter::CheckMargins(int imgwidth1, int imgheight1, int IcX, int IcY,
                           int ID)
{
    double LMarginMetricDeduct = 0.0;
    double RMarginMetricDeduct = 0.0;
    double TMarginMetricDeduct = 0.0;
    double BMarginMetricDeduct = 0.0;
    double OverallMarginDeduct = 0.0;
    LRMarginMetricTotalDeduct = 0.0;
    TBMarginMetricTotalDeduct = 0.0;

    LMargin = (double)(IcX - (ID / 2)) / (double)(ID / 2.0);
    RMargin = (double)(imgwidth1 - IcX - (ID / 2)) / (double)(ID / 2);
    TMargin = (double)(IcY - (ID / 2)) / (double)(ID / 2.0);
    BMargin = (double)(imgheight1 - IcY - (ID / 2)) / (double)(ID / 2);

    if (TMargin < 0.20) {
        // Deductions for Top Margin - 5/6 makes 1.0 Deduction at -1 LMargin
        TMarginMetricDeduct = (1.0 - (5.0 * TMargin)) / 6.0;
    }
    if (BMargin < 0.20) {
        BMarginMetricDeduct = (1.0 - (5.0 * BMargin)) / 6.0;
    }
    if (TMarginMetricDeduct < BMarginMetricDeduct) {
        TBMarginMetricTotalDeduct = BMarginMetricDeduct;
    }
    else {
        TBMarginMetricTotalDeduct = TMarginMetricDeduct;
    }
    if (TBMarginMetricTotalDeduct >= 1.0) {
        TBMarginMetricTotalDeduct = 1.0;
    }
    if (LMargin < 0.60) {
        // Deductions for Left Margin - 5/8 makes 1.0 Deduction at -1 LMargin
        LMarginMetricDeduct = (3.0 - (5.0 * LMargin)) / 8.0;
    }
    if (RMargin < 0.20) {
        RMarginMetricDeduct = (3.0 - (5.0 * RMargin)) / 8.0;
    }
    if (LMarginMetricDeduct < RMarginMetricDeduct) {
        LRMarginMetricTotalDeduct = RMarginMetricDeduct;
    }
    else {
        LRMarginMetricTotalDeduct = LMarginMetricDeduct;
    }
    if (LRMarginMetricTotalDeduct >= 1.0) {
        LRMarginMetricTotalDeduct = 1.0;
    }

    OverallMarginDeduct = TBMarginMetricTotalDeduct + LRMarginMetricTotalDeduct;

    if (OverallMarginDeduct > 1.0) {
        OverallMarginDeduct = 1.0;
    }

    OverallMargin = 1.0 - OverallMarginDeduct;
}

/**
 *
 * @param rawimg
 * @param imgwidth
 * @param imgheight
 * @param downsizescale
 * @param dsimage
 */
void MFilter::DownSize(unsigned char **rawimg, int imgwidth, int imgheight,
                       int downsizescale, unsigned char **dsimage)
{
    int indy;
    int indxx;
    int indyy;
    int origx;
    int origy;
    int dspoint;
    int averageboxsize = downsizescale * downsizescale;

    dswidth = imgwidth / downsizescale;
    dsheight = imgheight / downsizescale;

    for (int indx = 0; indx < dswidth; indx++) {
        origx = indx * downsizescale;
        for (indy = 0; indy < dsheight; indy++) {
            origy = indy * downsizescale;
            // for each pixel in the downsized image...
            // get the average intensity over the downsizescale*downsizescale
            // group of original pixels
            dspoint = 0;
            for (indxx = 0; indxx < downsizescale; indxx++) {
                for (indyy = 0; indyy < downsizescale; indyy++) {
                    dspoint += (int)rawimg[indxx + origx][indyy + origy];
                }
            }
            dsimage[indx][indy] = (unsigned char)(dspoint / averageboxsize);
        }
    }
}

/**
 *
 * @param imagein
 * @param dsimageinwidth
 * @param dsimageinheight
 */
void MFilter::EdgeMap(unsigned char **imagein, int dsimageinwidth,
                      int dsimageinheight)
{
    int indx;
    int jdx;
    int kdx;
    int imgpt1 = 0;
    int imgpt2 = 0;

    for (int indy = 4; indy < dsimageinheight - 5; indy++) {
        for (kdx = 0; kdx < 4; kdx++) {
            // assign 0 edge intensity values to a 4 pixel wide boundary on the
            // x axis note that this equates to a 16 pixel boundary on the
            // original image since it has been downsized by factor of 4
            Edgevalbin[kdx][indy] = 0;
            Edgevalbin[dsimageinwidth - 5 + kdx][indy] = 0;

            PupilEdgebin[kdx][indy] = 0;
            PupilEdgebin[dsimageinwidth - 5 + kdx][indy] = 0;

            EdgeRawbin[kdx][indy] = 0;
            EdgeRawbin[dsimageinwidth - 5 + kdx][indy] = 0;
        }
        for (indx = 4; indx < dsimageinwidth - 4; indx++) {
            Edgevalbin[indx][indy] = 0;
            PupilEdgebin[indx][indy] = 0;
            for (jdx = 0; jdx < 4;
                 jdx++) // Takes 4 pixels to either side of the trial pixel
            // Note that when this is done on a quarter downsized image - it is
            // a 16 pixel wide average difference
            {
                imgpt1 = (int)imagein[indx - jdx - 1][indy];
                imgpt2 = (int)imagein[indx + jdx][indy];
                // if (imgpt1 > 70) puppt1 = 70;
                // if (imgpt2 > 70) puppt2 = 70;
                if (imgpt1 > 200)
                    imgpt1 = 200; // Reduce the effects of specularities
                if (imgpt2 > 200)
                    imgpt2 = 200;
                Edgevalbin[indx][indy] += imgpt1 - imgpt2;
            }
            if (Edgevalbin[indx][indy] < 0) {
                Edgevalbin[indx][indy] *= -1;
            }
            if (Edgevalbin[indx][indy] > 140) {
                Edgevalbin[indx][indy] = 140;
            }
            // Edge differences greater than 140 are likely due to
            // specularities, not iris edges
            Edgevalbin[indx][indy] /= 3;
            EdgeRawbin[indx][indy] = (unsigned char)Edgevalbin[indx][indy];
        }
    }
}

/**
 *
 * @param dsbin
 * @param width
 * @param height
 */
void MFilter::GetVertEdges(unsigned char **dsbin, int width, int height)
{
    // Go horizontally across the downsized image and get average intensity
    // differences along an 8 pixel vertical edge
    // centered on the current horizontal point on each line
    int idx8;
    int jdy8;
    int kdy;
    int pointval;

    for (idx8 = 0; idx8 < width; idx8++) {
        // Assign zero vertical edge values in the first 4 vertical lines
        for (kdy = 0; kdy < 4; kdy++) {
            // assign 0 edge intensity values to a 4 pixel wide boundary on the
            // x axis note that this equates to a 16 pixel boundary on the
            // original image since it has been downsized by factor of 4
            VertEdge[idx8][kdy] = 0;
            VertEdge[idx8][height - kdy - 1] = 0;
        }
        for (jdy8 = 4; jdy8 < height - 8; jdy8++) {
            pointval = (int)dsbin[idx8][jdy8 - 2] + (int)dsbin[idx8][jdy8 - 1];
            pointval =
                pointval - (int)dsbin[idx8][jdy8 + 1] - (int)dsbin[idx8][jdy8];
            if (pointval < 0) {
                pointval = -pointval;
            }
            if (pointval > 255) {
                pointval = 255;
            }
            // Version 2.2.1 change
            // VertEdgebin[idx8][jdy8] = (unsigned char) pointval;
            VertEdge[idx8][jdy8] = (unsigned char)pointval;
            if (pointval > 20) {
                pointval = 20;
            }
            VertEdge[idx8][jdy8] = pointval;
        }
    }
}

/**
 *
 * @param Rawb
 * @param ImageWidth
 * @param ImageHeight
 * @param roughIxc
 * @param roughIyc
 * @param roughIdiam
 */
void MFilter::FindFineIris(unsigned char **Rawb, int ImageWidth,
                           int ImageHeight, int roughIxc, int roughIyc,
                           int roughIdiam)
{
    if (MathTablesDefined == false) {
        InitMathTables();
    }
    // set the bounds of the search area based upon rough locations
    int MinRfine = roughIdiam / 2 - 4;
    int MaxRfine = roughIdiam / 2 + 4;
    if (MaxRfine > MaxR) {
        MaxRfine = MaxR;
    }
    int MinIxcfine = roughIxc - 4;
    int MaxIxcfine = roughIxc + 4;
    int MinIycfine = roughIyc - 4;
    int MaxIycfine = roughIyc + 4;
    bestresponse = 0;
    leftresponse = 0;
    rightresponse = 0;
    int wf = 1;
    Point TanPts[33];

    int rightedgepointsused;
    int leftedgepointsused;
    int trialresponse;
    int leftedgeresponse;
    int rightedgeresponse;
    int srx;
    int suy;
    int slx;
    int sly;
    int idx;
    int jdx;
    int idy;
    int rctr;
    int jjdx;
    for (int rdx = MinRfine; rdx < MaxRfine; rdx++) {
        // calculate width and height of the test area for a given rdx, jdx
        // using real full resolution coordinates  int MaxRelX = rdx;  Duh... test
        // cannot include X coordinates within the trial radius  int MinRelY =
        // trigvals[LowerSegTrial, rdx+8, 1];  int MaxRelY =
        // trigvals[UpperSegTrial, rdx+8, 1];  remember that image y vals are
        // opposite to trig vals
        for (idx = MinIxcfine; idx < MaxIxcfine; idx++) {
            for (idy = MinIycfine; idy < MaxIycfine; idy++) {
                rightedgepointsused = 0;
                leftedgepointsused = 0;
                trialresponse = 0;
                leftedgeresponse = 0;
                rightedgeresponse = 0;
                // if (idx == RoughIriscx && idy == RoughIriscy && roughIdiam ==
                // rdx * 2)

                LowerSegTrial = 40; // 47 works  256 / 4 = 64
                UpperSegTrial = 30; // 35 works higher fails on droopy eyelids
                                    // cause failure lower fails on large pupils

                for (jdx = 0; jdx < LowerSegTrial; jdx += 2) // jdx++, jdx++)
                {
                    uppersegpointresponses[jdx] = 0;
                    // fix for the fact that the radial will extend an
                    // additional 12-16 pixels! Right extent
                    srx = idx + trigvals[jdx][rdx + 16][0];
                    // Upper extent - remember lower y is upper
                    suy = idy - trigvals[jdx][rdx + 16][1];
                    // Left extent
                    slx = idx - trigvals[jdx][rdx + 16][0];
                    // Lower extent
                    sly = idy + trigvals[jdx][rdx + 16][1];
                    wf = 1; // weighting factor to give more weight to pixels
                            // closest to the exact edge
                    if ((srx < ImageWidth) && (slx >= 0)) // Check to see if
                                                          // this is correct on
                                                          // the full size image
                    {
                        if ((suy >= 0) && (jdx < UpperSegTrial)) {
                            pointresponse = 0;
                            for (rctr = 0; rctr < 16; rctr++, rctr++) {
                                wf = 4 - (rctr / 4);
                                pointresponse +=
                                    wf *
                                    ((int)Rawb[idx +
                                               trigvals[jdx][rdx + rctr][0]]
                                              [idy -
                                               trigvals[jdx][rdx + rctr][1]] -
                                     (int)Rawb[idx +
                                               trigvals[jdx][rdx - rctr][0]]
                                              [idy -
                                               trigvals[jdx][rdx - rctr][1]]);
                            }
                            if (pointresponse < 0) {
                                pointresponse = -pointresponse;
                            }
                            if (pointresponse > maxpointresponse) {
                                pointresponse = maxpointresponse;
                            }
                            uppersegpointresponses[jdx] = pointresponse;
                            rightedgeresponse += pointresponse;
                            rightedgepointsused++;
                        }
                        if ((sly < ImageWidth) && (jdx < LowerSegTrial)) {
                            pointresponse = 0;
                            for (rctr = 0; rctr < 16; rctr++) {
                                wf = 4 - (rctr / 4);
                                pointresponse +=
                                    wf *
                                    ((int)Rawb[idx +
                                               trigvals[jdx][rdx + rctr][0]]
                                              [idy +
                                               trigvals[jdx][rdx + rctr][1]] -
                                     (int)Rawb[idx +
                                               trigvals[jdx][rdx - rctr][0]]
                                              [idy +
                                               trigvals[jdx][rdx - rctr][1]]);
                                // watch out for X values that go beyond the
                                // image
                            }
                            if (pointresponse < 0) {
                                pointresponse = -pointresponse;
                            }
                            if (pointresponse > maxpointresponse) {
                                pointresponse = maxpointresponse;
                            }
                            rightedgeresponse += pointresponse;
                            rightedgepointsused++;
                        }
                    }
                    if ((slx >= 0) && (srx <= ImageWidth)) {
                        if ((suy >= 0) && (jdx < UpperSegTrial)) {
                            pointresponse = 0;
                            for (rctr = 0; rctr < 16; rctr++) {
                                wf = 4 - (rctr / 4);
                                pointresponse +=
                                    wf *
                                    ((int)Rawb[idx -
                                               trigvals[jdx][rdx + rctr][0]]
                                              [idy -
                                               trigvals[jdx][rdx + rctr][1]] -
                                     (int)Rawb[idx -
                                               trigvals[jdx][rdx - rctr][0]]
                                              [idy -
                                               trigvals[jdx][rdx - rctr][1]]);
                            }
                            if (pointresponse < 0) {
                                pointresponse = -pointresponse;
                            }
                            if (pointresponse > maxpointresponse) {
                                pointresponse = maxpointresponse;
                            }
                            uppersegpointresponses[jdx] += pointresponse;
                            leftedgeresponse += pointresponse;
                            leftedgepointsused++;
                        }
                        if ((sly < ImageHeight - 1) && (jdx < LowerSegTrial)) {
                            // fill in logic to get the lower left response
                            pointresponse = 0;
                            for (rctr = 0; rctr < 16; rctr++) {
                                wf = 4 - (rctr / 4);
                                pointresponse +=
                                    wf *
                                    ((int)Rawb[idx -
                                               trigvals[jdx][rdx + rctr][0]]
                                              [idy +
                                               trigvals[jdx][rdx + rctr][1]] -
                                     (int)Rawb[idx -
                                               trigvals[jdx][rdx - rctr][0]]
                                              [idy +
                                               trigvals[jdx][rdx - rctr][1]]);
                            }
                            if (pointresponse < 0) {
                                pointresponse = -pointresponse;
                            }
                            if (pointresponse > maxpointresponse) {
                                pointresponse = maxpointresponse;
                            }
                            leftedgeresponse += pointresponse;
                            leftedgepointsused++;
                        }
                    }
                } // end of loop for angle segments to be added up
                trialresponse = leftedgeresponse + rightedgeresponse;

                if (trialresponse > bestresponse) {
                    bestresponse = trialresponse;

                    IrisCenterX = idx;
                    IrisCenterY = idy;
                    IrisRadius = rdx;

                    // now save the upperedgepointresponses for later use
                    for (jjdx = 0; jjdx < UpperSegTrial; jjdx++) {
                        maxuppersegresponses[jjdx] =
                            uppersegpointresponses[jjdx];
                    }
                }
            } // end of  y loop
        }     // end of x loop
    }         // end of R loop
    // Use the best diameter to calculate diagnostic metrics:
    // Get Iris-Sclera GrayScale Difference values over a -60 to 45 deg arc
    // along the diameter Use response value trends along the arc on left and
    // right sides to determine likely eyelid locations Fill in  ISGSVals[256];
    // // 0 is 0 degrees  // 64 is 90 deg.  //128 is 180 deg //192 is 270 deg.
    int IXC = IrisCenterX;
    int IYC = IrisCenterY;
    int IRAD = IrisRadius;

    memset(ISGSVals, 0, sizeof(int) * maxsegs);
    if (IrisRadius > MaxR - 16) {
        throw new std::out_of_range("Iris Radius exceeds max allowable");
    }
    if (IrisRadius < 17) {
        throw new std::out_of_range("Invalid Iris Diameter value");
    }

    int rctr2;
    int outerpointval;
    int innerpointval;
    for (int jdxx = 0; jdxx < maxsegs; jdxx++) {
        // in this array TanPts[33] is the centered on the iris edge arc
        // TanPts[16].X = IXC + trigvals[jdxx, IRAD, 0];
        // TanPts[16].Y = IXC + trigvals[jdxx, IRAD, 1];
        for (rctr = 0; rctr < 33; rctr++) {
            TanPts[rctr].X = IXC + trigvals[jdxx][IRAD + rctr - 16][0];
            TanPts[rctr].Y = IYC - trigvals[jdxx][IRAD + rctr - 16][1];
        }
        // if ((jdxx == 0) || (jdxx == 43) || (jdxx == 64) || (jdxx == 128) ||
        // (jdxx == 150))
        //{
        //    int checkpoint = 0;
        //}
        pointresponse = 0;
        int numradpts = 0;
        for (rctr2 = 1; rctr2 < 16; rctr2++) {
            if (TanPts[16 + rctr2].X > 0 && TanPts[16 + rctr2].X < ImageWidth &&
                TanPts[16 + rctr2].Y > 0 &&
                TanPts[16 + rctr2].Y < ImageHeight) {
                wf = 1;
                // weighting factor for distance from iris edge point going 16
                // pixels each direction wf = 4 - (rctr2 / 4);
                outerpointval =
                    (int)Rawb[TanPts[16 + rctr2].X][TanPts[16 + rctr2].Y];
                innerpointval =
                    (int)Rawb[TanPts[16 - rctr2].X][TanPts[16 - rctr2].Y];
                // Version 2.2.1 Change to limit effects of single differences)
                // pointresponse += wf * (outerpointval - innerpointval);
                int responseint = outerpointval - innerpointval;
                if (responseint > 75) {
                    responseint = 75;
                }
                pointresponse = pointresponse + responseint;
                numradpts++;
            }
            else {
                pointresponse = 0;
            }
        }
        ISGSVals[jdxx] = pointresponse;
    } // end of jdxx - intensity difference for each radial segment from 0 to
      // maxsegs at IXC, IYC, IRAD
    // Calculate adjusted Iris-Sclera GrayScale Difference using averages of
    // ISGSVals[128:138] and ISGSVals[246:0];
    int isgstotal = 0;
    // for (int jctr2 = 0; jctr2 < 10; jctr2++)
    //{
    //    isgstotal += ISGSVals[255 - jctr2];
    //    isgstotal += ISGSVals[128 + jctr2];
    //}
    // ISGSDiffMeanAvg = isgstotal / 300;   // remember that each point in
    // ISGSVals is 15 pixels on each side
    int ptsused = 0;
    for (int jctr2 = 0; jctr2 < 40; jctr2++) {
        // Each point in ISGSVals is 15 pixels on each side
        isgstotal = isgstotal + ISGSVals[255 - jctr2];
        ptsused++;
        isgstotal = isgstotal + ISGSVals[128 + jctr2];
        ptsused++;
    }
    for (int jctr3 = 0; jctr3 < 20; jctr3++) {
        // Each point in ISGSVals is 15 pixels on each side
        // Add response on 20 upper segments
        isgstotal = isgstotal + ISGSVals[jctr3];
        ptsused++;
        isgstotal = isgstotal + ISGSVals[127 - jctr3];
        ptsused++;
    }
    ISGSDiffMeanAvg =
        isgstotal /
        (ptsused *
         15); // Remember that each point in ISGSVals is 15 pixels on each side
}

/**
 *
 * @param framebytes
 * @param width
 * @param height
 */
void MFilter::FindIris(const unsigned char *framebytes, int width, int height)
{
    int i;
    unsigned char **MyTwoDArray;

    MyTwoDArray = new unsigned char *[width];
    for (i = 0; i < width; i++) {
        MyTwoDArray[i] = new unsigned char[height];
    }

    RawByteFrameToTwoDimArray(framebytes, MyTwoDArray, width, height);
    FindIris(MyTwoDArray, width, height);

    for (i = 0; i < width; i++) {
        delete[] MyTwoDArray[i];
    }
    delete[] MyTwoDArray;
}

/**
 *
 * @param RawbytesIn
 * @param ImageWidth
 * @param ImageHeight
 */
void MFilter::FindIris(unsigned char **RawbytesIn, int ImageWidth,
                       int ImageHeight)
{
    int i;
    int j;
    int dsfactor = 4;
    int dswidth = ImageWidth / dsfactor;
    int dsheight = ImageHeight / dsfactor;

    //    MBitmap nbmp = new MBitmap();
    // MBMP.RawBMP nrbmp = new MBMP.RawBMP();
    // nbmp.WriteRawBytesToBMPFile(DSRawbin,ImageWidth/4,
    // ImageHeight/4,@"C:\temp\test.bmp"
    for (i = 0; i < ImageWidth; i++) {
        for (j = 0; j < ImageHeight; j++) {
            Rawbin[i][j] = RawbytesIn[i][j];
        }
    }

    DownSize(RawbytesIn, ImageWidth, ImageHeight, dsfactor, DSRawbin);
    // nbmp.WriteRawBytesToBMPFile(DSRawbin,ImageWidth/4,
    // ImageHeight/4,@"C:\temp\DSImage.bmp");
    EdgeMap(DSRawbin, dswidth, dsheight);
    GetVertEdges(DSRawbin, dswidth, dsheight);
    FindIrisCenter(Edgevalbin, Edgevalbin, dswidth, dsheight, Rawbin);
}

/**
 *
 * @param Ledgevals
 * @param Redgevals
 * @param inimagewidth
 * @param inimageheight
 * @param Rawbytes
 */
void MFilter::FindIrisCenter(int **Ledgevals, int **Redgevals, int inimagewidth,
                             int inimageheight, unsigned char **Rawbytes)
{
    if (MathTablesDefined == false) {
        InitMathTables();
    }
    // Version 2.2.1 change initialize darkiriscenter to zero
    int darkiriscenter = 0;
    bestdarkiriscenter = 0;

    Point iriscenterpt[224];
    // assume that Edgevalbin from above is used
    // remember that image coordinants are from top left to bottom right
    // MinRtrial = 35;  //Reminder that since this is on a half-sized image...
    // correlates with 150 iris diameter min MaxRtrial = 100;  // Max trial
    // radius for 320x240 image... Max trial diameter is 400 on full-sized
    // 640x480 image LowerSegTrial = 24;   is the number of 360/256 degree
    // segments below the horizontal axis of the iris UpperSegTrial = 12;   is
    // the number of segments above the horizontal axis of the iris

    // Version 2.2.1 Change ... LowerSegTrial and UpperSegTrial are globals...
    //                          Reset these values due to conflict with Fine
    //                          Iris finding routine.
    LowerSegTrial = 47;
    UpperSegTrial = 40;

    int MaxSegTrial = UpperSegTrial;
    if (LowerSegTrial > MaxSegTrial) {
        MaxSegTrial = LowerSegTrial;
    }

    bestresponse = 0;
    RoughIrisRadius = 0;
    RoughIriscx = 0;
    RoughIriscy = 0;
    IrisCenterX = 0;
    IrisCenterY = 0;
    IrisRadius = 0;
    // set up the array of circular points close to the center to check for dark
    // areas. set up the upper right quadrant and copy to the others
    int ctr = 0;
    int rdxx;
    int jdxx;
    for (rdxx = 2; rdxx < 17;
         rdxx +=
         2) // makes points up to 16 pixel (downsized) radius (8 samples)
    // equivalent pupil diameter is 16*2 * 4 = 128
    {
        for (jdxx = 0; jdxx < 61;
             jdxx +=
             10) // takes 10 angle segment increments up to 60 ( 7 samples )
        {
            // creates 56 points in upper right trig quadrant (in relative
            // coordinates, not image coords)
            iriscenterpt[ctr].X = trigvals[jdxx][rdxx][0];
            iriscenterpt[ctr].Y = trigvals[jdxx][rdxx][1];
            // creates additional 56 points in lower right, upper left and lower
            // left quadrants
            iriscenterpt[ctr + 56].X = trigvals[jdxx][rdxx][0];
            iriscenterpt[ctr + 56].Y = -trigvals[jdxx][rdxx][1];
            iriscenterpt[ctr + 112].X = -trigvals[jdxx][rdxx][0];
            iriscenterpt[ctr + 112].Y = trigvals[jdxx][rdxx][1];
            iriscenterpt[ctr + 168].X = -trigvals[jdxx][rdxx][0];
            iriscenterpt[ctr + 168].Y = -trigvals[jdxx][rdxx][1];
            ctr++;
        }
    }
    maxIrisFindRightEdgeResponse = 0;
    maxIrisFindLeftEdgeResponse = 0;

    int rdx;
    int idx;
    int idy;
    int jdx;
    int MaxRelX;
    int MinRelY;
    int MaxRelY;
    int rightedgepointsused;
    int leftedgepointsused;
    int trialresponse;
    int leftedgeresponse;
    int rightedgeresponse;
    int vertedgeresponse;
    int ptdarkval;
    int trialx;
    int trialy;
    int ictr5;
    int srx;
    int suy;
    int slx;
    int sly;
    for (rdx = MinRtrial; rdx < MaxRtrial; rdx++, rdx++) {
        // calculate width and height of the test area for a given rdx,jdx
        MaxRelX =
            rdx; // test cannot include X coordinants within the trial radius
        MinRelY = trigvals[LowerSegTrial][rdx][1];
        MaxRelY = trigvals[UpperSegTrial][rdx][1];
        // remember that these are opposite image coords
        for (idx = inimagewidth * 3 / 10; idx < inimagewidth * 7 / 10;
             idx++, idx++) {
            for (idy = inimageheight * 3 / 10; idy < inimageheight * 7 / 10;
                 idy++, idy++) {
                rightedgepointsused = 0;
                leftedgepointsused = 0;
                trialresponse = 0;
                leftedgeresponse = 0;
                rightedgeresponse = 0;
                vertedgeresponse = 0;
                darkiriscenter = 0;
                ptdarkval = 0;
                trialx = 0;
                trialy = 0;
                for (ictr5 = 0; ictr5 < 224; ictr5++) {
                    // calculate "Dark Center response" for the current idx,idy
                    // note that the grid radius is always 5,10,15,20,25
                    trialx = idx + iriscenterpt[ictr5].X;
                    trialy = idy - iriscenterpt[ictr5].Y;
                    if (trialx > 0 && trialx < inimagewidth &&
                        trialy < inimageheight && trialy > 0) {
                        ptdarkval = (int)DSRawbin[trialx][trialy];
                        // ptdarkval = Rawbin[trialx, trialy];
                        ptdarkval = (55 - ptdarkval);
                        // ptdarkval = (130 - ptdarkval);  gives too much credit
                        // in low contrast images to dark gray iris
                        if (ptdarkval < 0) {
                            ptdarkval = 0;
                        }
                        // ptdarkval = (ptdarkval ) ;
                        darkiriscenter += ptdarkval;
                    }
                }
                for (jdx = 0; jdx < MaxSegTrial; jdx++) {
                    // right extent
                    srx = idx + trigvals[jdx][rdx][0];
                    // upper extent - remember lower y is upper
                    suy = idy - trigvals[jdx][rdx][1];
                    // left extent
                    slx = idx - trigvals[jdx][rdx][0];
                    // lower extent
                    sly = idy + trigvals[jdx][rdx][1];

                    if ((srx + 4 < inimagewidth - 1) && (srx - 2 >= 0)) {
                        if ((suy >= 4) && (jdx < UpperSegTrial)) {
                            rightedgeresponse += Redgevals[srx][suy];
                            vertedgeresponse += VertEdge[srx][suy];
                            rightedgeresponse += Edgevalbin[srx - 1][suy];
                            vertedgeresponse += VertEdge[srx - 1][suy];
                            rightedgeresponse += Edgevalbin[srx + 1][suy];
                            vertedgeresponse += VertEdge[srx + 1][suy];
                            rightedgepointsused++;
                        }
                        if ((sly < inimageheight - 3) &&
                            (jdx < LowerSegTrial)) {
                            rightedgeresponse += Edgevalbin[srx][sly];
                            vertedgeresponse += VertEdge[srx][sly];
                            rightedgeresponse += Edgevalbin[srx - 1][sly];
                            vertedgeresponse += VertEdge[srx - 1][sly];
                            rightedgeresponse += Edgevalbin[srx + 1][sly];
                            vertedgeresponse += VertEdge[srx + 1][sly];
                            rightedgepointsused++;
                        }
                    }

                    // Left Edge response
                    if ((slx - 2 >= 0) && (srx + 4 < inimagewidth - 1)) {
                        if ((suy >= 4) && (jdx < UpperSegTrial)) {
                            leftedgeresponse += Edgevalbin[slx][suy];
                            vertedgeresponse += VertEdge[slx][suy];
                            leftedgeresponse += Edgevalbin[slx - 1][suy];
                            vertedgeresponse += VertEdge[slx - 1][suy];
                            leftedgeresponse += Edgevalbin[slx + 1][suy];
                            vertedgeresponse += VertEdge[slx + 1][suy];
                            leftedgepointsused++;
                        }
                        if ((sly < inimageheight - 1) &&
                            (jdx < LowerSegTrial)) {
                            // lower left edge response
                            leftedgeresponse += Edgevalbin[slx][sly];
                            vertedgeresponse += VertEdge[slx][sly];
                            leftedgeresponse += Edgevalbin[slx - 1][sly];
                            vertedgeresponse += VertEdge[slx - 1][sly];
                            leftedgeresponse += Edgevalbin[slx + 1][sly];
                            vertedgeresponse += VertEdge[slx + 1][sly];
                            leftedgepointsused++;
                        }
                    }
                } // end of loop for angle segments to be added up
                // weighting factor for darkiriscenter
                darkiriscenter =
                    darkiriscenter * 105 /
                    10; // 100 used in test - works with "dark threshold" of 55
                vertedgeresponse = vertedgeresponse * 3 / 10;
                trialresponse = leftedgeresponse * 15 / 10 +
                                rightedgeresponse * 15 / 10 + vertedgeresponse +
                                darkiriscenter;

                if (trialresponse > bestresponse) {
                    bestresponse = trialresponse;
                    RoughIriscx = idx * 4;
                    RoughIriscy = idy * 4;
                    RoughIrisRadius = rdx * 4;
                    IrisCenterX = idx * 4;
                    IrisCenterY = idy * 4;
                    IrisRadius = rdx * 4;
                    maxIrisFindRightEdgeResponse = rightedgeresponse;
                    maxIrisFindLeftEdgeResponse = leftedgeresponse;
                    maxVertEdgeResponse = vertedgeresponse;
                    bestdarkiriscenter = darkiriscenter;
                }
            } // end of  y loop
        }     // end of x loop
    }         // end of R loop
    // check for strong vertical edges along the iris-sclera boundary as an
    // indicator of occlusions
    FindFineIris(Rawbin, inimagewidth * 4, inimageheight * 4, IrisCenterX,
                 IrisCenterY, IrisRadius * 2);
}

/**
 *
 * @param Edgev
 * @param inimagewidth
 * @param inimageheight
 * @param iriscx
 * @param iriscy
 * @param irisd
 */
void MFilter::FindPupilCenter(int **Edgev, int inimagewidth, int inimageheight,
                              int iriscx, int iriscy, int irisd)
{
    try {
        // Use same routine as Find Iris but restrict search to +/-0.1Idiam
        // inside iris region and limit radius range to 0.8 irisd. This means
        // that pupil center is restricted to within 40% of Iris diameter from
        // Iris center. Also restrict trial iris radius so that it does not get
        // within 8 pixels from iris edge
        int irisreduced = irisd / 10;
        int roitop = iriscy - irisreduced;
        int roileft = iriscx - irisreduced;
        int roibottom = iriscy + irisreduced;
        int roiright = iriscx + irisreduced;
        int PupIcdist = 0;
        int ipgsresponse = 0;
        MaxPupilEdgeResponse = 0;
        edgeresponse = 0;
        int pwf = 0;
        Point pupilcenterpt[224];

        int pctr = 0;
        int rdxx;
        int jdxx;
        for (rdxx = 3; rdxx < 25;
             rdxx +=
             3) // makes points up to 45 radius = 100 pixel diameter (8 samples)
        // this is in downsized (by factor of 4) dimensions
        {
            // takes 10 angle segment increments up to 60 ( 7 samples )
            for (jdxx = 0; jdxx < 61; jdxx += 10)
            {
                // creates 56 points in upper right trig quadrant (in relative
                // coordinates, not image coords)
                pupilcenterpt[pctr].X = trigvals[jdxx][rdxx][0];
                pupilcenterpt[pctr].Y = trigvals[jdxx][rdxx][1];
                // creates additional 56 points in lower right, upper left and
                // lower left quadrants
                pupilcenterpt[pctr + 56].X = trigvals[jdxx][rdxx][0];
                pupilcenterpt[pctr + 56].Y = -trigvals[jdxx][rdxx][1];
                pupilcenterpt[pctr + 112].X = -trigvals[jdxx][rdxx][0];
                pupilcenterpt[pctr + 112].Y = trigvals[jdxx][rdxx][1];
                pupilcenterpt[pctr + 168].X = -trigvals[jdxx][rdxx][0];
                pupilcenterpt[pctr + 168].Y = -trigvals[jdxx][rdxx][1];
                pctr++;
            }
        }
        int ipgspointsused = 0;
        // temp test

        PUpperSegTrial = 60;
        PLowerSegTrial = 63;
        // make sure ROI is inside the image
        if (roitop < 0) {
            roitop = 0;
        }
        if (roileft < 0) {
            roileft = 0;
        }
        if (roiright > inimagewidth - 1) {
            roiright = inimagewidth - 1;
        }
        if (roibottom > inimageheight - 1) {
            roibottom = inimageheight - 1;
        }
        int PMaxSegTrial = PUpperSegTrial;
        if (PLowerSegTrial > PMaxSegTrial) {
            PMaxSegTrial = PLowerSegTrial;
        }
        // int darkpupilresponse = 0;
        int bestdarkpupilresponse = 0;
        int darkpupilval = 0;
        int bestresponse = 0;
        // int pptresponse = 0;
        MaxPupilEdgeResponse = 0;
        edgeresponse = 0;
        pupilCx = 0;
        pupilCy = 0;
        // for each trial radius rdx  and for angle range jdx... find the best
        // total response for the image

        MaxPrtrial = IrisRadius / 4 * 8 / 10; // Use 80% of IrisRadius
        int rightedgepointsused = 0;
        int leftedgepointsused = 0;
        int rdx;
        int idx;
        int idy;
        int jdx;
        int ptdarkval;
        int trialx;
        int trialy;
        int pctr5;
        int trialresponse;
        int leftedgeresponse;
        int rightedgeresponse;
        int srx;
        int suy;
        int slx;
        int sly;
        for (rdx = MinPrtrial; rdx < MaxPrtrial; rdx++) {
            // calculate width and height of the test area for a given rdx,jdx
            // int MaxRelX = rdx; //Test cannot include X coordinants within the
            // trial radius  int MinRelY = trigvals[LowerSegTrial, rdx, 1];  int
            // MaxRelY = trigvals[UpperSegTrial, rdx, 1];  //remember that these
            // are opposite image coords

            for (idx = roileft; idx < roiright; idx++) {
                for (idy = roitop; idy < roibottom; idy++) {
                    darkpupilval = 0;
                    ptdarkval = 0;
                    trialx = 0;
                    trialy = 0;
                    for (pctr5 = 0; pctr5 < 224; pctr5++) {
                        // calculate "Dark Center response" for the current idx,
                        // idy note that the grid radius is always 5,10,15,20,25
                        trialx = idx + pupilcenterpt[pctr5].X;
                        trialy = idy - pupilcenterpt[pctr5].Y;
                        if (trialx > 0 && trialx < inimagewidth &&
                            trialy < inimageheight && trialy > 0) {
                            ptdarkval = (int)DSRawbin[trialx][trialy];
                            ptdarkval = 130 - ptdarkval;
                            if (ptdarkval < 0) {
                                ptdarkval = 0;
                            }
                            ptdarkval /= 5;
                            darkpupilval = darkpupilval + ptdarkval;
                        }
                    }
                    // darkpupilresponse = darkpupilresponse / 20;
                    ipgsresponse = 0;
                    ipgspointsused = 0;
                    rightedgepointsused = 0;
                    leftedgepointsused = 0;
                    trialresponse = 0;
                    leftedgeresponse = 0;
                    rightedgeresponse = 0;

                    // for (int jdx = 0; jdx < MaxSegTrial;  jdx++, jdx++ )
                    for (jdx = 0; jdx < PMaxSegTrial; jdx++, jdx++) {
                        srx = idx + trigvals[jdx][rdx][0]; // Right extent
                        suy = idy - trigvals[jdx][rdx][1]; // Upper extent -
                                                           // remember lower y
                                                           // is upper
                        slx = idx - trigvals[jdx][rdx][0]; // Left extent
                        sly = idy + trigvals[jdx][rdx][1]; // Lower extent

                        if ((srx + 4 < inimagewidth - 1) && (srx - 4 >= 0)) {
                            if ((suy >= 0) && (jdx < PUpperSegTrial)) {
                                // rightedgeresponse = rightedgeresponse + 1*
                                // Redgevals[srx, suy];
                                rightedgeresponse += 2 * Edgev[srx][suy];
                                // rightedgeresponse = rightedgeresponse + 1 *
                                // Redgevals[srx - 2, suy];
                                rightedgeresponse += Edgev[srx - 1][suy];
                                rightedgeresponse += Edgev[srx + 1][suy];
                                // rightedgeresponse = rightedgeresponse + 1 *
                                // Redgevals[srx + 2, suy];
                                rightedgepointsused++;
                            }
                            if ((sly < inimageheight - 1) &&
                                (jdx < PLowerSegTrial)) {
                                rightedgeresponse += 2 * Edgev[srx][sly];
                                // rightedgeresponse = rightedgeresponse + 1 *
                                // Redgevals[srx - 2, sly];
                                rightedgeresponse += Edgev[srx - 1][sly];
                                rightedgeresponse += Edgev[srx + 1][sly];
                                // rightedgeresponse = rightedgeresponse + 1 *
                                // Redgevals[srx + 2, sly];
                                rightedgepointsused++;
                            }
                        }
                        if ((slx - 4 >= 0) && (srx + 4 <= inimagewidth - 1)) {
                            if ((suy >= 0) && (jdx < PUpperSegTrial)) {
                                leftedgeresponse += 2 * Edgev[slx][suy];
                                // leftedgeresponse = leftedgeresponse + 1 *
                                // Ledgevals[slx - 2, suy];
                                leftedgeresponse += Edgev[slx - 1][suy];
                                leftedgeresponse += Edgev[slx + 1][suy];
                                // leftedgeresponse = leftedgeresponse + 1 *
                                // Ledgevals[slx + 2, suy];
                                leftedgepointsused++;
                            }
                            if ((sly < inimageheight - 1) &&
                                (jdx < PLowerSegTrial)) {
                                leftedgeresponse += 2 * Edgev[slx][sly];
                                // leftedgeresponse = leftedgeresponse + 1 *
                                // Ledgevals[slx - 2, sly];
                                leftedgeresponse += Edgev[slx - 1][sly];
                                leftedgeresponse += Edgev[slx + 1][sly];
                                // leftedgeresponse = leftedgeresponse + 1 *
                                // Ledgevals[slx + 2, sly];
                                leftedgepointsused++;
                            }
                        }
                    } // end of loop for angle segments to be added up
                    if (leftedgepointsused + rightedgepointsused) {
                        edgeresponse = (leftedgeresponse + rightedgeresponse) *
                                       100 /
                                       (leftedgepointsused + rightedgepointsused);
                    } else {
                        edgeresponse = 0;
                    }
                    darkpupilval *= 2; // 4
                    trialresponse = edgeresponse + darkpupilval;
                    if (trialresponse > bestresponse) {
                        bestresponse = trialresponse;
                        MaxPupilEdgeResponse = edgeresponse;
                        bestdarkpupilresponse = darkpupilval;
                        pupilCx = (idx - 1) * 4;
                        pupilCy = (idy)*4;
                        pupilRad = (rdx)*4;
                        pupilD = rdx * 8;
                    }
                    {
                        PupIcdist = (int)sqrt(
                            (double)((iriscy - idy) * (iriscy - idy) +
                                     (iriscx - idx) * (iriscx - idx)));
                        if ((PupIcdist + rdx) > (irisd / 2) * 65 / 100 ||
                            (idy + rdx) > (iriscy + (irisd / 2) * 65 / 100) ||
                            (idy - rdx) < (idy - (iriscy / 2) * 65 / 100))
                            // combination of trial radius rdx and trial pupilCx
                            // (idx), pupilCy (idy) is within 85% of the edge of
                            // the Iris so stop the loop
                            idy = roibottom;
                    }
                } // end of  y loop

                // PupIcdist = (int)Math.Sqrt((double)((iriscy - idy) * (iriscy
                // - idy) + (iriscx - idx) * (iriscx - idx)));
                if ((idx + rdx) > (iriscx + (irisd / 2) * 65 / 100) ||
                    (iriscx - rdx) < (idx - (irisd / 2) * 65 / 100)) {
                    // combination of trial radius rdx and trial pupilCx (idx),
                    // pupilCy (idy) is within 85% of the edge of the Iris so
                    // stop the loop
                    idx = roiright;
                }
            } // end of x loop
        }     // end of R loop  ROUGH PUPIL RADIUS and Center point Done

        RoughPcx = pupilCx;
        RoughPcy = pupilCy;
        RoughPrad = pupilRad;
        RoughPdiam = pupilRad * 2;

        // calculate Pupil Quality Metrics
        ipgsresponse = 0;
        ipgspointsused = 0;
        int idxxx = pupilCx / 4;
        int idyyy = pupilCy / 4;
        int rdrrr = pupilRad / 4;
        int jdxxx;
        int jdyyy;
        for (jdxxx = -20; jdxxx < 20; jdxxx++) {
            for (jdyyy = -5; jdyyy < 5; jdyyy++) {
                if ((idxxx + jdxxx - rdrrr) > 0 &&
                    ((idxxx + jdxxx + rdrrr) < (inimagewidth - 1)) &&
                    (idyyy + jdyyy) > 0 &&
                    ((idyyy + jdyyy) < (inimageheight - 1))) {
                    ipgsresponse += Edgev[idxxx - rdrrr + jdxxx][idyyy + jdyyy];
                    ipgspointsused++;
                    ipgsresponse += Edgev[idxxx + rdrrr + jdxxx][idyyy + jdyyy];
                    ipgspointsused++;
                }
            }
        }
        IrisPupilGSDiff = (float)ipgsresponse / (float)ipgspointsused;
        PupilIrisDiameterRatio = ((float)pupilRad / (float)IrisRadius);
    }
    catch (std::exception &ex) {
        pupilCx = IrisCenterX;
        pupilCy = IrisCenterY;
        pupilRad = IrisRadius;
        pupilD = IrisRadius * 2;
        IrisPupilGSDiff = 0;
        PupilIrisDiameterRatio = (float)0.8;
        std::cerr << "exception caught: " << ex.what() << std::endl;
    }
}

/**
 *
 * @param imagein
 * @param imwidth
 * @param imheight
 * @param imicx
 * @param imicy
 * @param imirad
 * @param impcx
 * @param impcy
 * @param imprad
 */
void MFilter::FindOcclusions(unsigned char **imagein, int imwidth, int imheight,
                             int imicx, int imicy, int imirad, int impcx,
                             int impcy, int imprad)
{
    int irelx = 0;
    int irelxx = 0;
    int pupradsquared = imprad * imprad;
    int xtest = 0;
    int ytest = 0;
    int iyupper = 0;
    int iylower = 0;
    int cctr;
    int dctr;
    int ictr;
    int dx;
    int dy = 0;
    int xchordpt = 0;
    int edgelim = 80;
    int speclim = 225;
    int specdefval = 200;
    int ptedge;
    int refptedge;
    int ptval;

    float zval = 5.8f;

    int upperrightgray[100];
    int upperleftgray[100];
    int lowerleftgray[100];
    int lowerrightgray[100];
    int urgpts[100];
    int ulgpts[100];
    int lrgpts[100];
    int llgpts[100];
    int uloccpts[100];
    int uroccpts[100];
    int lloccpts[100];
    int lroccpts[100];

    Point ulgstart[100];
    Point urgstart[100];
    Point llgstart[100];
    Point lrgstart[100];
    Point ulgstop[100];
    Point urgstop[100];
    Point llgstop[100];
    Point lrgstop[100];

    memset(upperrightgray, 0, sizeof(int) * 100);
    memset(upperleftgray, 0, sizeof(int) * 100);
    memset(lowerleftgray, 0, sizeof(int) * 100);
    memset(lowerrightgray, 0, sizeof(int) * 100);
    memset(urgpts, 0, sizeof(int) * 100);
    memset(ulgpts, 0, sizeof(int) * 100);
    memset(lrgpts, 0, sizeof(int) * 100);
    memset(llgpts, 0, sizeof(int) * 100);
    memset(uroccpts, 0, sizeof(int) * 100);
    memset(uloccpts, 0, sizeof(int) * 100);
    memset(lroccpts, 0, sizeof(int) * 100);
    memset(lloccpts, 0, sizeof(int) * 100);
    memset(ulgstart, 0, sizeof(Point) * 100);
    memset(urgstart, 0, sizeof(Point) * 100);
    memset(llgstart, 0, sizeof(Point) * 100);
    memset(lrgstart, 0, sizeof(Point) * 100);
    memset(ulgstop, 0, sizeof(Point) * 100);
    memset(urgstop, 0, sizeof(Point) * 100);
    memset(llgstop, 0, sizeof(Point) * 100);
    memset(lrgstop, 0, sizeof(Point) * 100);

    for (cctr = 0; cctr < imwidth; cctr++) {
        for (dctr = 0; dctr < imheight; dctr++) {
            MaskedRawbin[cctr][dctr] = Rawbin[cctr][dctr];
        }
    }
    for (ictr = 0; ictr < 100; ictr++) {
        // dy is trigonometric vertical distance of the test line from the
        // centerline of iris (remember to subtract when used in image space) dy
        // is multiplied by 2 since each loop increases the vertical distance by
        // 2/100 (1/50) of the top half of the iris
        dy = ictr;
        irelxx = (int)floor(0.5 + sqrt((double)(10000 - dy * dy)));
        // corresponding maximum DX from center for circle with 100 radius
        // irelxx = (int)Math.Round(Math.Sqrt((double)(10000 - dy * dy))) ;

        // test case when ictr is 25 (50% height -> 30 deg.), irelxrightent =
        // (86.7 ) * 120 / 100 = 103  -> (dx,dy) =
        // (+103, +60) from center point of 120 rad circle
        // irisnormchords[25,0] is (-103,50)

        for (dx = 0; dx < irelxx - 8; dx++) // added -8  to avoid including
                                            // actual iris edge points as
                                            // "occlusions"
        {
            // test case  imagewidth = 640 imageheight = 480  imicx = 320  imicy
            // = 240  imirad=125  impcx = 330 impcy = 250   imprad = 40  test
            // case ictr = 25
            iyupper = imicy - imirad * dy / 100; //
            // iyupper = 177 = 240- 63(i.e.. the top scan line of the iris is at
            // iris center y - iris radius )
            iylower = imicy + imirad * dy / 100;
            // test case:  iylower = 303  = 240 + 63  at ictr = 25 (50% from
            // center to top edge)
            xchordpt =
                imicx +
                dx * imirad /
                    100; //    iyupper, iylower, irelx are actual locations
            // get right side iris intensities

            if (xchordpt < imwidth) // checking to make sure the x value is
                                    // within the right margin of the image
            {
                xtest =
                    xchordpt -
                    impcx; // xtest is the x distance from the pupil center x.
                if ((iyupper >= 0)) // checking that top of the iris does not
                                    // extend beyond the upper image boundary
                {
                    // check to see if the point is inside the pupil boundary

                    ytest = iyupper - impcy; // iyupper is the current y
                                             // coordinate value for the chord

                    // will not count anything within +/- 32 pixels of the
                    // assumed pupil radius outline
                    if ((xtest * xtest) + (ytest * ytest) >
                        pupradsquared + 20 * 20) {
                        if (xchordpt > 0 && xchordpt < imwidth && iyupper > 0 &&
                            iyupper < imheight - 1) {
                            ptval = (int)imagein[xchordpt][iyupper];
                            if (ptval > speclim)
                                ptval = specdefval;
                            upperrightgray[ictr] = upperrightgray[ictr] + ptval;
                            urgpts[ictr]++;
                            if (urgstart[ictr].X == 0) // record the leftmost
                                                       // pixel that is included
                                                       // in this scan line
                            {
                                urgstart[ictr].X = xchordpt;
                                urgstart[ictr].Y = iyupper;
                            }
                            // update the end pixel position...  the last one
                            // recorded will be used
                            urgstop[ictr].X = xchordpt;
                            urgstop[ictr].Y = iyupper;
                            // ValF16Edge  ValF16EdgeBytes
                            refptedge =
                                (int)VertEdgebin[xchordpt / 4][iyupper / 4] +
                                (int)EdgeRawbin[xchordpt / 4][iyupper / 4];
                            ptedge =
                                (int)VertEdgebin[xchordpt / 4][iyupper / 4] +
                                ValF16Edge[xchordpt][iyupper];
                            if ((int)imagein[xchordpt][iyupper] > speclim ||
                                ptedge > edgelim) {
                                MaskedRawbin[xchordpt][iyupper] = 0xFF;
                                // point is likely a specularity or is on an
                                // "edge" that indicates an eyelid
                                uroccpts[ictr]++;
                            }
                        }
                    }
                }
                if (iylower <= imheight - 1) {
                    ytest = iylower - impcy;
                    if ((xtest * xtest) + (ytest * ytest) >
                        pupradsquared + 20 * 20) {
                        if (xchordpt > 0 && xchordpt < imwidth && iylower > 0 &&
                            iylower < imheight) {
                            ptval = (int)imagein[xchordpt][iylower];
                            if (ptval > speclim)
                                ptval = specdefval;
                            lowerrightgray[ictr] = lowerrightgray[ictr] + ptval;
                            lrgpts[ictr]++;
                            if (lrgstart[ictr].X == 0) // record the leftmost
                                                       // pixel that is included
                                                       // in this scan line
                            {
                                lrgstart[ictr].X = xchordpt;
                                lrgstart[ictr].Y = iylower;
                            }
                            // update the end pixel position...  the last one
                            // recorded will be used
                            lrgstop[ictr].X = xchordpt;
                            lrgstop[ictr].Y = iylower;
                            ptedge =
                                (int)VertEdgebin[xchordpt / 4][iylower / 4] +
                                (int)EdgeRawbin[xchordpt / 4][iylower / 4];
                            if (imagein[xchordpt][iylower] > speclim ||
                                ptedge > edgelim) {
                                MaskedRawbin[xchordpt][iylower] = 0xFF;
                                // point is likely a specularity or is on an
                                // "edge" that indicates an eyelid
                                lroccpts[ictr]++;
                            }
                        }
                    }
                }
            }
        } // finished with Right Side

        for (dx = 0; dx < irelxx - 8; dx++) // added -8 to avoid including
                                            // actual iris edge points as
                                            // "occlusions"
        {
            // test case  imagewidth = 640 imageheight = 480  imicx = 320  imicy
            // = 240  imirad=125  impcx = 330 impcy = 250   imprad = 40  test
            // case ictr = 25
            iyupper = imicy - imirad * dy / 100;
            // iyupper = 177 = 240- 63(i.e.. the top scan line of the iris is at
            // iris center y - iris radius )
            iylower = imicy + imirad * dy / 100;
            // test case:  iylower = 303  = 240 + 63  at ictr = 25 (50% from
            // center to top edge)
            xchordpt =
                imicx -
                dx * imirad /
                    100; //    iyupper, iylower, irelx are actual locations
            // get right side iris intensities

            if (xchordpt < imwidth - 1) // checking to make sure the x value is
                                        // within the right margin of the image
            {
                xtest = xchordpt - impcx;
                // checking that top of the iris does not extend beyond the
                // upper image boundary (y=0)
                if ((iyupper >= 0)) {
                    // check to see if the point is inside the pupil boundary

                    ytest = iyupper - impcy;

                    if ((xtest * xtest) + (ytest * ytest) >
                        pupradsquared + 20 * 20) {
                        if (xchordpt > 0 && xchordpt < imwidth && iyupper > 0 &&
                            iyupper < imheight) {
                            ptval = (int)imagein[xchordpt][iyupper];
                            if (ptval > speclim)
                                ptval = specdefval;
                            upperleftgray[ictr] = upperleftgray[ictr] + ptval;
                            ulgpts[ictr]++;
                            if (ulgstop[ictr].X == 0) // record the leftmost
                                                      // pixel that is included
                                                      // in this scan line
                            {
                                ulgstop[ictr].X = xchordpt;
                                ulgstop[ictr].Y = iyupper;
                            }
                            // update the end pixel position...  the last one
                            // recorded will be used
                            ulgstart[ictr].X = xchordpt;
                            ulgstart[ictr].Y = iyupper;
                            ptedge =
                                (int)VertEdgebin[xchordpt / 4][iyupper / 4] +
                                (int)EdgeRawbin[xchordpt / 4][iyupper / 4];
                            if (imagein[xchordpt][iyupper] > speclim ||
                                ptedge > edgelim) {
                                MaskedRawbin[xchordpt][iyupper] = 0xFF;
                                // point is likely a specularity or is on an
                                // "edge" that indicates an eyelid
                                uloccpts[ictr]++;
                            }
                        }
                    }
                }
            }
            if (iylower <= imheight - 1) {
                ytest = iylower - impcy;
                if ((xtest * xtest) + (ytest * ytest) >
                    pupradsquared + 20 * 20) {
                    if (xchordpt > 0 && xchordpt < imwidth && iylower > 0 &&
                        iylower < imheight) {
                        ptval = (int)imagein[xchordpt][iylower];
                        if (ptval > speclim)
                            ptval = specdefval;
                        lowerleftgray[ictr] = lowerleftgray[ictr] + ptval;
                        llgpts[ictr]++;
                        if (llgstop[ictr].X == 0) // record the leftmost pixel
                                                  // that is included in this
                                                  // scan line
                        {
                            llgstop[ictr].X = xchordpt;
                            llgstop[ictr].Y = iylower;
                        }
                        // update the end pixel position...  the last one
                        // recorded will be used
                        llgstart[ictr].X = xchordpt;
                        llgstart[ictr].Y = iylower;
                        ptedge = (int)VertEdgebin[xchordpt / 4][iylower / 4] +
                                 (int)EdgeRawbin[xchordpt / 4][iylower / 4];
                        if (imagein[xchordpt][iylower] > speclim ||
                            ptedge > edgelim) {
                            MaskedRawbin[xchordpt][iylower] = 0xFF;
                            // point is likely a specularity or is on an "edge"
                            // that indicates an eyelid
                            lloccpts[ictr]++;
                        }
                    }
                }
            }
        }
    }
    for (int kctr = 0; kctr < 100; kctr++) {
        if (urgpts[kctr] > 0) {
            upperrightgray[kctr] /= urgpts[kctr];
        }
        if (lrgpts[kctr] > 0) {
            lowerrightgray[kctr] /= lrgpts[kctr];
        }
        if (ulgpts[kctr] > 0) {
            upperleftgray[kctr] /= ulgpts[kctr];
        }
        if (llgpts[kctr] > 0) {
            lowerleftgray[kctr] /= llgpts[kctr];
        }
    }
    int lowerrightreference = 0;
    int lowerleftreference = 0;
    int mctr3;
    for (int mctr = 0; mctr < 50; mctr++) {
        lowerrightreference += lowerrightgray[mctr];
        lowerleftreference += lowerleftgray[mctr];
    }
    int lowerrightrefavg = lowerrightreference / 50;
    int lowerleftrefavg = lowerleftreference / 50;
    int lowerrightvar = 0;
    int lowerleftvar = 0;
    for (int mctr2 = 5; mctr2 < 50; mctr2++) {
        lowerrightvar += (lowerrightgray[mctr2] - lowerrightrefavg) *
                         (lowerrightgray[mctr2] - lowerrightrefavg);
        lowerleftvar += (lowerleftgray[mctr2] - lowerleftrefavg) *
                        (lowerleftgray[mctr2] - lowerleftrefavg);
    }
    int combinedlowref = 0;
    for (mctr3 = 5; mctr3 < 45; mctr3++) {
        combinedlowref += lowerrightgray[mctr3] + lowerleftgray[mctr3];
    }
    int combinedvar = 0;
    int combinedlowavg = combinedlowref / ((40) * 2);
    for (int mctr4 = 5; mctr4 < 45; mctr4++) {
        combinedvar += (lowerrightgray[mctr4] - combinedlowavg) *
                       (lowerrightgray[mctr4] - combinedlowavg);
        combinedvar += (lowerleftgray[mctr4] - combinedlowavg) *
                       (lowerleftgray[mctr4] - combinedlowavg);
    }
    int stddevcombined =
        (int)floor(0.5 + sqrt((double)(combinedvar) / (double)(90)));
    // int stddevright = (int)(Math.Round(Math.Sqrt((double)lowerrightvar /
    // (double)40)));  int stddevleft =
    // (int)(Math.Round(Math.Sqrt((double)lowerleftvar / (double)40)));
    int leftavailpixels = 0;
    int rightavailpixels = 0;
    int lefttotalpixels = 0;
    int righttotalpixels = 0;
    int lrblock3 = 0;
    int urblock3 = 0;
    int llblock3 = 0;
    int ulblock3 = 0;
    int actr;
    for (mctr3 = 0; mctr3 < 100; mctr3++) {
        // find how much is occluded
        if (lrblock3 < 3) {
            // lower right quadrant
            if (abs(lowerrightgray[mctr3] - lowerrightrefavg) >
                (int)((float)stddevcombined * zval)) {
                lrblock3++;
                // area is occluded
                righttotalpixels += lrgpts[mctr3];
                // mark it in the MaskedRawbin

                for (actr = lrgstart[mctr3].X; actr < lrgstop[mctr3].X;
                     actr++) {
                    MaskedRawbin[actr][lrgstart[mctr3].Y] = 0;
                }
            }
            else {
                lrblock3 = 0; // if there are less than 3 non-occluded lines in
                              // a row... reset the
                righttotalpixels += lrgpts[mctr3];
                rightavailpixels += lrgpts[mctr3] - lroccpts[mctr3];
            }
        }
        else {
            // area follows three occluded lines in a row..  Mark all remaining
            // lines as occluded.
            righttotalpixels += lrgpts[mctr3];
            // mark it in the MaskedRawbin

            for (actr = lrgstart[mctr3].X; actr < lrgstop[mctr3].X; actr++) {
                MaskedRawbin[actr][lrgstart[mctr3].Y] = 0;
            }
        }

        // upper right
        if (urblock3 < 3) {
            if (abs(upperrightgray[mctr3] - lowerrightrefavg) >
                (int)((float)stddevcombined * zval)) {
                // area is occluded
                urblock3++;
                righttotalpixels += urgpts[mctr3];
                // mark it in the MaskedRawbin

                for (actr = urgstart[mctr3].X; actr < urgstop[mctr3].X;
                     actr++) {
                    MaskedRawbin[actr][urgstart[mctr3].Y] = 0x00;
                }
            }
            else {
                righttotalpixels += urgpts[mctr3];
                rightavailpixels += urgpts[mctr3] - uroccpts[mctr3];
                urblock3 = 0; // reset the "three in a row" counter
            }
        }
        else {
            // area follows three occluded lines in a row..  mark all remaining
            // lines as occluded.
            righttotalpixels += urgpts[mctr3];
            // mark it in the MaskedRawbin

            for (actr = urgstart[mctr3].X; actr < urgstop[mctr3].X; actr++) {
                MaskedRawbin[actr][urgstart[mctr3].Y] = 0;
            }
        }
        // lower left quadrant
        if (llblock3 < 3) {
            if (abs(lowerleftgray[mctr3] - lowerleftrefavg) >
                (int)((float)stddevcombined * zval)) {
                // area is occluded
                llblock3++;
                lefttotalpixels += llgpts[mctr3];
                // mark it in the MaskedRawbin

                for (actr = llgstart[mctr3].X; actr < llgstop[mctr3].X;
                     actr++) {
                    MaskedRawbin[actr][llgstart[mctr3].Y] = 0;
                }
            }
            else {
                lefttotalpixels += llgpts[mctr3];
                leftavailpixels += llgpts[mctr3] - lloccpts[mctr3];
                llblock3 = 0;
            }
        }
        else {
            lefttotalpixels += llgpts[mctr3];
            // mark it in the MaskedRawbin

            for (actr = llgstart[mctr3].X; actr < llgstop[mctr3].X; actr++) {
                MaskedRawbin[actr][llgstart[mctr3].Y] = 0;
            }
        }

        // upper left
        if (ulblock3 < 3) {
            if (abs(upperleftgray[mctr3] - lowerleftrefavg) >
                (int)((float)stddevcombined * zval)) {
                // area is occluded
                ulblock3++;
                lefttotalpixels += ulgpts[mctr3];
                // mark it in the MaskedRawbin

                for (actr = ulgstart[mctr3].X; actr < ulgstop[mctr3].X;
                     actr++) {
                    MaskedRawbin[actr][ulgstart[mctr3].Y] = 0x00;
                }
            }
            else {
                lefttotalpixels += ulgpts[mctr3];
                leftavailpixels += ulgpts[mctr3] - uloccpts[mctr3];
                ulblock3 = 0;
            }
        }
        else {
            // area is occluded
            lefttotalpixels += ulgpts[mctr3];
            // mark it in the MaskedRawbin

            for (actr = ulgstart[mctr3].X; actr < ulgstop[mctr3].X; actr++) {
                MaskedRawbin[actr][ulgstart[mctr3].Y] = 0;
            }
        }
    }
    int leftavailpercent = 0;
    if (lefttotalpixels > 0) {
        leftavailpercent *= 100 / lefttotalpixels;
    }
    int rightavailpercent = 0;
    if (righttotalpixels > 0) {
        rightavailpercent *= 100 / righttotalpixels;
    }
    int combinedavailpercent = 0;
    if (lefttotalpixels + righttotalpixels > 0) {
        combinedavailpercent = 100 * (leftavailpixels + rightavailpixels) /
                               (lefttotalpixels + righttotalpixels);
    }
    if (combinedavailpercent > 100) {
        combinedavailpercent = 100;
    }
    if (combinedavailpercent < 0) {
        combinedavailpercent = 0;
    }
    UsableIrisAreaPercent = combinedavailpercent;
}

/**
 *
 */
void MFilter::InitMathTables()
{
    int anglesegnum;
    int radiuslength;
    int halfmaxsegs = maxsegs / 2;
    int seg90 = (int)((double)maxsegs / (double)4);
    int seg180 = (int)((double)maxsegs / (double)2);
    double dsin;
    double dcos;
    double singlesegment = (double)((double)6.28318530718 / (double)maxsegs);

    for (anglesegnum = 0; anglesegnum < seg90 + 1; anglesegnum++) {
        // calculate sine and cosines
        dsin = sin(singlesegment * (double)(anglesegnum));
        dcos = ::cos(singlesegment * (double)(anglesegnum));
        for (radiuslength = 0; radiuslength < MaxR + 1; radiuslength++) {
            // get x axis orinals in image array for a given angle and distance
            // from center
            if (dcos >= 0) {
                trigvals[anglesegnum][radiuslength][0] =
                    (int)floor(dcos * (double)radiuslength + 0.5);
            }
            else {
                trigvals[anglesegnum][radiuslength][0] =
                    (int)ceil(dcos * (double)radiuslength - 0.5);
            }
            trigvals[halfmaxsegs - anglesegnum][radiuslength][0] =
                -trigvals[anglesegnum][radiuslength][0] - 1;
            trigvals[halfmaxsegs + anglesegnum][radiuslength][0] =
                -trigvals[anglesegnum][radiuslength][0] - 1;
            trigvals[maxsegs - anglesegnum][radiuslength][0] =
                trigvals[anglesegnum][radiuslength][0];

            // now get y axis ordinal in image array for a given angle and
            // distance from center
            if (dsin >= 0) {
                trigvals[anglesegnum][radiuslength][1] =
                    (int)floor(dsin * (double)radiuslength + 0.5);
            }
            else {
                trigvals[anglesegnum][radiuslength][1] =
                    (int)ceil(dsin * (double)radiuslength - 0.5);
            }
            trigvals[halfmaxsegs - anglesegnum][radiuslength][1] =
                trigvals[anglesegnum][radiuslength][1];
            trigvals[halfmaxsegs + anglesegnum][radiuslength][1] =
                -trigvals[anglesegnum][radiuslength][1] - 1;
            trigvals[maxsegs - anglesegnum][radiuslength][1] =
                -trigvals[anglesegnum][radiuslength][1];
        } // radius
    }
    // special cases just to keep from referencing out-of-bounds values at MaxR
    // since defined origin at lower left corner of (0,0)
    trigvals[0][MaxR][0] = MaxR - 1;
    trigvals[seg90][MaxR][1] = MaxR - 1;
    trigvals[halfmaxsegs][MaxR][0] = -MaxR;
    trigvals[maxsegs - seg90][MaxR][1] = -MaxR;
    trigvals[maxsegs][MaxR][0] = MaxR - 1;
    MathTablesDefined = true;
}

/**
 *
 */
void MFilter::FillEdgeTable()
{
    for (int ctr = 0; ctr < 100; ctr++) {
        //  each ctr++ is finds coordinants for the upper right edge of a circle
        //  with radius 100
        irispts[ctr].Y = ctr;
        irispts[ctr].X = (int)floor(
            0.5 + sqrt((double)10000 - (double)pow((double)(ctr), 2)));
    }
    irispts[100].X = 0;
    irispts[100].Y = 100;
    // irispts[100] is top dead center
}

/**
 *
 * @param contrast_in
 * @param defocus_in
 * @param ISGSMean_in
 * @param Margin_in
 * @param IrisD_in
 * @param IrisVis_in
 * @param ipgsdiff_in
 * @param pupilirisratio_in
 * @return
 */
int MFilter::CalcOverallQuality(int contrast_in, int defocus_in,
                                int ISGSMean_in, float Margin_in, int IrisD_in,
                                float IrisVis_in, float ipgsdiff_in,
                                float pupilirisratio_in)
{
    OverallQuality = 0;
    // Normalize Contrast value
    NormalizeContrast(contrast_in);
    // Normalize Defocus score
    NormalizeDefocus(defocus_in);
    // Normalize ISGS score
    NormalizeISGS(ISGSMean_in);
    // Normalize IrisDiamter
    NIrisID = NormalizeIrisDiameter(IrisD_in);
    // Normalize IPGSDifference
    NIPGSDiff = NormalizeIPGSDiff(ipgsdiff_in);
    // Normalize IPDiamRatio
    NPupilIrisRatio = NormalizePupilIrisRatio(pupilirisratio_in);
    // Normalize IrisVisibility
    NIrisVis = NormalizeIrisVisibility(IrisVis_in);
    // Calculate CombinedQuality Score
    OverallQuality =
        GetQuality(NContrast, NDefocus, NIrisID, NISGSMean, NIrisVis, Margin_in,
                   NPupilIrisRatio, NIPGSDiff);
    // Finally  Scale the CombinedQuality output metric
    return OverallQuality;
}

/**
 *
 * @param contrast_in
 * @return
 */
float MFilter::NormalizeContrast(int contrast_in)
{
    NContrast = 0;
    if (contrast_in < ContrastLowLimit) {
        NContrast = 0;
    }
    if (contrast_in >= ContrastLowLimit && contrast_in <= ContrastUpperLimit) {
        NContrast = (float)(((float)(contrast_in)-ContrastLowLimit) /
                            (ContrastUpperLimit - ContrastLowLimit));
    }
    if (contrast_in >= ContrastUpperLimit) {
        NContrast = 1;
    }
    return NContrast;
}

/**
 *
 * @param defocus_in
 * @return
 */
float MFilter::NormalizeDefocus(int defocus_in)
{
    NDefocus = 0;
    if (defocus_in < DefocusLowLimit) {
        NDefocus = 0;
    }
    if (defocus_in >= DefocusLowLimit && defocus_in <= DefocusUpperLimit) {
        NDefocus = (float)(((float)(defocus_in)-DefocusLowLimit) /
                           (DefocusUpperLimit - DefocusLowLimit));
    }
    if (defocus_in >= DefocusUpperLimit) {
        NDefocus = 1;
    }
    return NDefocus;
}

/**
 *
 * @param ISGSMean_in
 * @return
 */
float MFilter::NormalizeISGS(int ISGSMean_in)
{
    NISGSMean = 0;
    // if (ISGSMean_in < ISGSMeanLowLimit) NISGSMean = 0.8f;
    if (ISGSMean_in < ISGSMeanLowLimit) {
        NISGSMean = 0.0; // Revision 2.2 on April 10,2012
    }
    if (ISGSMean_in >= ISGSMeanLowLimit && ISGSMean_in <= ISGSMeanUpperLimit) {
        NISGSMean =
            (float)(0.8 + (ISGSMean_in - ISGSMeanLowLimit) /
                              (ISGSMeanUpperLimit - ISGSMeanLowLimit) * 0.2);
    }
    if (ISGSMean_in >= ISGSMeanUpperLimit) {
        NISGSMean = 1;
    }
    return NISGSMean;
}

/**
 *
 * @param IrisD_in
 * @return
 */
float MFilter::NormalizeIrisDiameter(int IrisD_in)
{
    NIrisID = 0;
    if (IrisD_in < IDLowLimit) {
        NIrisID = 0;
    }
    if (IrisD_in >= IDLowLimit && IrisD_in <= IDMedLimit1) {
        NIrisID = (IrisD_in - IDLowLimit) / (IDMedLimit1 - IDLowLimit);
    }
    if (IrisD_in >= IDMedLimit1 && IrisD_in <= IDMedLimit2) {
        NIrisID = 1.0;
    }
    if (IrisD_in >= IDMedLimit2 && IrisD_in <= IDUpperLimit) {
        NIrisID = (float)(1.0 - ((float)IrisD_in - IDMedLimit2) /
                                    (IDUpperLimit - IDMedLimit2));
    }
    if (IrisD_in >= IDUpperLimit) {
        NIrisID = 0;
    }
    return NIrisID;
}

/**
 *
 * @param ipgsdiff_in
 * @return
 */
float MFilter::NormalizeIPGSDiff(float ipgsdiff_in)
{
    float NIPGSDiffValueAtLowLimit = (float)0.85;
    NIPGSDiff = 0;
    if (ipgsdiff_in < IPGSDiffLowLimit) {
        NIPGSDiff = (float)0.85;
    }
    if (ipgsdiff_in >= IPGSDiffLowLimit && ipgsdiff_in <= IPGSDiffUpperLimit) {
        NIPGSDiff = (float)(NIPGSDiffValueAtLowLimit +
                            (ipgsdiff_in - IPGSDiffLowLimit) /
                                (IPGSDiffUpperLimit - IPGSDiffLowLimit) *
                                (1.0 - NIPGSDiffValueAtLowLimit));
    }
    if (ipgsdiff_in >= IPGSDiffUpperLimit) {
        NIPGSDiff = 1;
    }
    return NIPGSDiff;
}

/**
 *
 * @param pupilirisratio_in
 * @return
 */
float MFilter::NormalizePupilIrisRatio(float pupilirisratio_in)
{
    NPupilIrisRatio = 0.0;
    if (pupilirisratio_in < PupilIrisRatioLowLimit) {
        NPupilIrisRatio = 0.0;
    }
    if (pupilirisratio_in >= PupilIrisRatioLowLimit &&
        pupilirisratio_in <= PupilIrisRatioHighLimit) {
        NPupilIrisRatio = 1.0;
    }
    if (pupilirisratio_in > PupilIrisRatioHighLimit) {
        NPupilIrisRatio = 0.0;
    }
    return NPupilIrisRatio;
}

/**
 *
 * @param IrisVis_in
 * @return
 */
float MFilter::NormalizeIrisVisibility(float IrisVis_in)
{
    NIrisVis = 0;
    if (IrisVis_in < IrisVisLowLimit) {
        NIrisVis = 0;
    }
    if (IrisVis_in >= IrisVisLowLimit && IrisVis_in <= IrisVisUpperLimit) {
        NIrisVis =
            (float)(0.7 + (IrisVis_in - IrisVisLowLimit) /
                              (IrisVisUpperLimit - IrisVisLowLimit) * 0.3);
    }
    if (IrisVis_in >= IrisVisUpperLimit) {
        NIrisVis = 1;
    }
    return NIrisVis;
}

/**
 *
 * @param NContrast
 * @param NDefocus
 * @param NIrisID
 * @param NISGSMean
 * @param NIrisVis
 * @param NMargin_in
 * @param NPupilIrisRatio_in
 * @param NIPGSDiff
 * @return
 */
int MFilter::GetQuality(float NContrast, float NDefocus, float NIrisID,
                        float NISGSMean, float NIrisVis, float NMargin_in,
                        float NPupilIrisRatio_in, float NIPGSDiff)
{
    CombinedQuality =
        (float)(NContrast * NContrastFactor * NDefocus * NDefocusFactor *
                NIrisID * NISGSMean * NISGSFactor * NIrisVis * NIrisVisFactor *
                NMargin_in * NPupilIrisRatio_in * NPupilIrisRatioFactor *
                NIPGSDiff * NIPGSDiffFactor);

    if (CombinedQuality > CombinedQualityLowLimit) {
        OverallQuality = 0;
    }
    if (CombinedQuality >= CombinedQualityLowLimit &&
        CombinedQuality <= CombinedQualityUpperLimit) {
        OverallQuality =
            (int)(80.0 * (CombinedQuality - CombinedQualityLowLimit) /
                  (CombinedQualityUpperLimit - CombinedQualityLowLimit));
    }
    if (CombinedQuality > CombinedQualityUpperLimit) {
        OverallQuality =
            (int)(80.0 + 20.0 * (CombinedQuality - CombinedQualityUpperLimit) /
                             (1.0 - CombinedQualityUpperLimit));
    }
    return OverallQuality;
}

/**
 *
 * @param framebytes
 * @param width
 * @param height
 * @return
 */
int MFilter::fastQuality(const unsigned char *framebytes, int width, int height)
{
    if (width > maxwidth || height > maxheight) {
        std::cerr << "ERROR: Maximum image dimensions exceeded. Images should "
                     "be smaller than "
                  << maxwidth << "x" << maxheight << "." << std::endl;
        return ERROR_IMAGESIZE_MAX_EXCEEDED;
    }
    RawByteFrameToTwoDimArray(framebytes, Rawbin, width, height);
    int DefResult = Defocus(Rawbin, width, height);
    int ConResult = Contrast(Rawbin, width, height);
    NDefocus = NormalizeDefocus(DefResult);
    NContrast = NormalizeContrast(ConResult);
    return GetQuality(NContrast, NDefocus, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
}

/**
 *
 * @param rawbytes
 * @param width
 * @param height
 * @return
 */
int MFilter::FastQuality(unsigned char **rawbytes, int width, int height)
{
    int DefResult = Defocus(Rawbin, width, height);
    int ConResult = Contrast(Rawbin, width, height);
    NDefocus = NormalizeDefocus(DefResult);
    NContrast = NormalizeContrast(ConResult);
    return GetQuality(NContrast, NDefocus, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
}

/**
 *
 * @param framein
 * @param width
 * @param height
 * @return
 */
int MFilter::getQualityFromImageFrame(const unsigned char *framein, int width,
                                      int height)
{
    FindIris(framein, width, height);
    Defocus(Rawbin, width, height);
    Contrast(Rawbin, width, height);
    // IrisScleraGSDiff(Rawbin, width, height, IrisCenterX, IrisCenterY,
    // IrisRadius * 2);
    CheckMargins(width, height, IrisCenterX, IrisCenterY, IrisRadius * 2);
    FindPupilCenter(Edgevalbin, width / 4, height / 4, IrisCenterX / 4,
                    IrisCenterY / 4, IrisRadius * 2 / 4);
    FindPupilCenter(Edgevalbin, width / 4, height / 4, IrisCenterX / 4,
                    IrisCenterY / 4, IrisRadius * 2 / 4);
    FindOcclusions(Rawbin, width, height, IrisCenterX, IrisCenterY, IrisRadius,
                   pupilCx, pupilCy, pupilRad);
    float irisvis = (float)((float)UsableIrisAreaPercent / 10.0);
    return CalcOverallQuality(ContrastScore, DefocusScore, (int)ISGSDiffMeanAvg,
                              (float)OverallMargin, IrisRadius * 2, irisvis,
                              IrisPupilGSDiff, PupilIrisDiameterRatio);
}

/**
 *
 * @param rawin
 * @param width
 * @param height
 * @return
 */
int MFilter::GetQualityFromImage(unsigned char **rawin, int width, int height)
{
    FindIris(rawin, width, height);
    Defocus(Rawbin, width, height);
    Contrast(Rawbin, width, height);
    // IrisScleraGSDiff(Rawbin, width, height, IrisCenterX, IrisCenterY,
    // IrisRadius * 2);
    CheckMargins(width, height, IrisCenterX, IrisCenterY, IrisRadius * 2);
    FindPupilCenter(Edgevalbin, width / 4, height / 4, IrisCenterX / 4,
                    IrisCenterY / 4, IrisRadius * 2 / 4);
    FindPupilCenter(Edgevalbin, width / 4, height / 4, IrisCenterX / 4,
                    IrisCenterY / 4, IrisRadius * 2 / 4);
    FindOcclusions(Rawbin, width, height, IrisCenterX, IrisCenterY, IrisRadius,
                   pupilCx, pupilCy, pupilRad);
    float irisvis = (float)((float)UsableIrisAreaPercent / 10.0);
    return CalcOverallQuality(ContrastScore, DefocusScore, (int)ISGSDiffMeanAvg,
                              (float)OverallMargin, IrisRadius * 2, irisvis,
                              IrisPupilGSDiff, PupilIrisDiameterRatio);
}

/**
 *
 */
void MFilter::initialize()
{
    if (MathTablesDefined == false) {
        InitMathTables();
        FillEdgeTable();
    }
}
