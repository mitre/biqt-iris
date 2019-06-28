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
//        FILE: ImageOps.h
// DESCRIPTION: Class that provides open source Defocus
//              and Overall Contrast quality metrics for
//              .bmp iris images.
//        DATE: 13 June 2012
//
//---------------------------------------------------------------------

#ifndef IMAGE_OPS_H
#define IMAGE_OPS_H

#ifdef _WIN32
#define IMAGEOPS_EXPORT __declspec(dllexport)
#else
#define IMAGEOPS_EXPORT
#endif

#define FOCUS_KERNEL_DIM 8
#define HALF_FOCUS_KERNEL_DIM 4
#define ERROR_NO_ROI_SELECTED -1001
#define ERROR_IMAGESIZE_MAX_EXCEEDED -1002
#define ERROR_ZERO_IMAGESIZE -1003
#define ERROR_NON_8BPP_INDEXED_BITMAP_FILE -1004
#define ERROR_LEFT_SCLERA_MARGIN -1101
#define ERROR_RIGHT_SCLERA_MARGIN -1102
#define ERROR_UNSPECIFIED_ERROR -9999
#define RIGHT_IRIS_GS_MARGIN_VALID true
#define LEFT_IRIS_GS_MARGIN_VALID true
#define BOTTOM_IRIS_GS_MARGIN_VALID true
#define TOP_IRIS_GS_MARGIN_VALID true

struct Point {
    int X;
    int Y;
};

// Class for basic image operations
// Built for  Defocus Metric and Overall Contrast Metric
class IMAGEOPS_EXPORT MFilter {
  public:
    MFilter();
    ~MFilter();
    void initialize();
    int fastQuality(const unsigned char *framebytes, int width, int height);
    int getQualityFromImageFrame(const unsigned char *framein, int width,
                                 int height);
    inline int getOverallQuality();
    inline int getDefocusScore();
    inline int getContrastScore();
    inline int getIrisRadius();
    inline int getUsableIrisAreaPercent();
    inline int getIrisCenterX();
    inline int getIrisCenterY();
    inline float getNDefocus();
    inline float getNContrast();
    inline float getNISGSMean();
    inline float getIrisPupilGSDiff();
    inline float getNIPGSDiff();
    inline float getNIrisID();
    inline float getNIrisVis();
    inline double getISGSDiffMeanAvg();

  private:
    void RawByteFrameToTwoDimArray(const unsigned char *singlearray,
                                   unsigned char **twodimarray, int width,
                                   int height);
    int Defocus(unsigned char **rawimg, int imgwidth, int imgheight);
    int Defocus(unsigned char *framebytes, int width, int height);
    int Contrast(unsigned char **rawimg2, int imgwidth2, int imgheight2);
    int Contrast(unsigned char *framebytes, int imgwidth, int imgheight);
    void CheckMargins(int imgwidth1, int imgheight1, int IcX, int IcY, int ID);
    void DownSize(unsigned char **rawimg, int imgwidth, int imgheight,
                  int downsizescale, unsigned char **dsimage);
    void EdgeMap(unsigned char **imagein, int dsimageinwidth,
                 int dsimageinheight);
    void GetVertEdges(unsigned char **dsbin, int width, int height);
    void FindFineIris(unsigned char **Rawb, int ImageWidth, int ImageHeight,
                      int roughIxc, int roughIyc, int roughIdiam);
    void FindIris(const unsigned char *framebytes, int width, int height);
    void FindIris(unsigned char **RawbytesIn, int ImageWidth, int ImageHeight);
    void FindIrisCenter(int **Ledgevals, int **Redgevals, int inimagewidth,
                        int inimageheight, unsigned char **Rawbytes);
    void FindPupilCenter(int **Edgev, int inimagewidth, int inimageheight,
                         int iriscx, int iriscy, int irisd);
    void FindOcclusions(unsigned char **imagein, int imwidth, int imheight,
                        int imicx, int imicy, int imrad, int impcx, int impcy,
                        int imprad);
    void InitMathTables();
    void FillEdgeTable();
    int CalcOverallQuality(int contrast_in, int defocus_in, int ISGSMean_in,
                           float Margin_in, int IrisD_in, float IrisVis_in,
                           float ipgsdiff_in, float pupilirisratio_in);
    int GetQuality(float NContrast, float NDefocus, float NIrisID,
                   float NISGSMean, float NIrisVis, float NMargin_in,
                   float NPupilIrisRatio_in, float NIPGSDiff);
    int FastQuality(unsigned char **rawbytes, int width, int height);
    int GetQualityFromImage(unsigned char **rawin, int width, int height);
    float NormalizeContrast(int contrast_in);
    float NormalizeDefocus(int defocus_in);
    float NormalizeISGS(int ISGSMean_in);
    float NormalizeIrisDiameter(int IrisD_in);
    float NormalizeIPGSDiff(float ipgsdiff_in);
    float NormalizePupilIrisRatio(float pupilirisratio_in);
    float NormalizeIrisVisibility(float IrisVis_in);

    static const int maxwidth = 640;  // RMD will want to change/increase
    static const int maxheight = 480; // these
    static const int FocusKernel[FOCUS_KERNEL_DIM][FOCUS_KERNEL_DIM];
    static const int HalfFocusKernel[HALF_FOCUS_KERNEL_DIM]
                                    [HALF_FOCUS_KERNEL_DIM];

    bool MathTablesDefined;

    int DefocusScore;
    int DefocusScoreRef;
    float NormDefocusScore;

    unsigned char *Framebin;
    unsigned char **Rawbin;
    int **ValF16Edge;
    unsigned char **ValF16EdgeBytes;
    unsigned char **TransposedRawbin;
    unsigned char **MaskedRawbin;
    int **Rawvals;
    unsigned char **DSRawbin;
    int **DSvalbin;
    unsigned char **NegDSRawbin;
    int **NegDSvalbin;
    int **Edgevalbin;
    int **PupilEdgebin;
    int **TransposedEdgebin;
    unsigned char **EdgeRawbin;
    int **VertEdge;
    unsigned char **VertEdgebin;

    int dswidth;
    int dsheight;
    static const int maxsegs = 256; // use a lower number for quicker execution
    static const int MaxR = 250;    // use a lower number for quicker execution
    int IrisCenterX;
    int IrisCenterY;
    int IrisRadius;
    int RoughIrisRadius;
    int RoughIriscx;
    int RoughIriscy;
    int RoughPcx;
    int RoughPcy;
    int RoughPrad;
    int RoughPdiam;
    int bestresponse;
    int rightresponse;
    int leftresponse;
    int maxIrisFindLeftEdgeResponse;
    int bestdarkiriscenter;
    int maxIrisFindRightEdgeResponse;
    int MaxPupilEdgeResponse;
    int edgeresponse;
    int maxVertEdgeResponse;
    int leftedgeIrisCenterx;
    int leftedgeIrisCentery;
    int leftedgeIrisDiameter;
    int rightedgeIrisCenterx;
    int rightedgeIrisCentery;
    int rightedgeIrisDiameter;
    int UsableIrisAreaPercent;
    static int PolarSegMax;

    /** The angle in the range 0-255. 0 is 0 degrees. 64 is 90 deg. 128 is 180
     * deg 192 is 270 deg */
    int *ISGSVals;
    int *PupilVals;

    float ContrastLowLimit;
    float ContrastUpperLimit;
    float NContrastValueAtLowLimit;
    float NContrastValueAtUpperLimit;
    float NContrastFactor;
    float DefocusLowLimit;
    float DefocusUpperLimit;
    float NDefocusValueAtLowLimit;
    float NDefocusValueAtUpperLimit;
    float NDefocusFactor;
    float IDLowLimit;
    float IDMedLimit1;
    float IDMedLimit2;
    float IDUpperLimit;
    float NIDiamFactor;
    float ISGSMeanLowLimit;
    float ISGSMeanUpperLimit;
    float NISGSValueAtLowLimit;
    float NISGSValueAtUpperLimit;
    float NISGSFactor;
    float IrisVisLowLimit;
    float IrisVisUpperLimit;
    float IPGSDiffLowLimit;
    float IPGSDiffUpperLimit;
    float NIPGSDiffFactor;
    float NPupilIrisRatio;
    float PupilIrisRatioLowLimit;
    float PupilIrisRatioHighLimit;
    float NPupilIrisRatioFactor;
    float NIrisVisFactor;
    float CombinedQuality;
    float NContrast;
    float NDefocus;
    float NIrisID;
    float NISGSMean;
    float NIrisVis;
    float NIPGSDiff;
    float CombinedQualityLowLimit;
    float CombinedQualityUpperLimit;

    int OverallQuality;
    int ***trigvals;
    /**
     * holds x coords. for 1% vertical steps up a 100 radius circle for upper
     * right quadrant ex.:  irispts[50].Y = 100*sin(30deg)
     */
    Point *irispts;
    int *cos;

    /**
     * sets up a table of trig values where [t, , ] is the anglesegment theta in
     * 2pi radians/maxsegs increments [ ,r, ] is the radius distance [ , ,0] is
     * r*sin(theta) [ , ,1] is r*cos(theta)
     */

    /**
     * Iris Finder Settings
     *  27 is good ..25 is 100 radius - 200 diameter 22 is 88 radius - 176 diam
     *  42 is good for test data (max diameter = 336) 400 diameter - 400/8 = 50
     *  60 is good Max trial radius for 160x120 image...60 is 240 radius - 480
     * diameter
     */
    int MinRtrial;
    int MaxRtrial;
    /** 47 works  256 / 4 = 64  43 is 60 degrees   32 is 45 deg. */
    static int LowerSegTrial;
    /** 35 works higher fails on droopy eyelids cause failure   lower fails on
     * large pupils */
    static int UpperSegTrial;
    /** Upper max segment MUST ALWAYS BE LESS THAN number of lowersegments tried
     */
    int *uppersegpointresponses;
    int *maxuppersegresponses;
    /** Pupil finder settings */
    static int PLowerSegTrial;
    /** Pupil finder settings */
    static int PUpperSegTrial;

    /**
     * 6 x4 x 2 = 48 pixel diameter    8 x 4*2 = 64 pixel diameter
     * This value will get changed by the pupil finder based upon iris diameter
     * found 22 * 8 is 176 pixel max pupil diamter
     */
    int MinPrtrial;
    int MaxPrtrial;
    int pupilCx;
    int pupilCy;
    int pupilRad;
    int pupilD;
    int ***radials;

    static double dMaxR;
    static int seg90;
    static double singlesegment;

    int ContrastScore;
    int AverageImageIntensity;

    double LMargin;
    double RMargin;
    double TMargin;
    double BMargin;
    double LRMarginMetricTotalDeduct;
    double TBMarginMetricTotalDeduct;
    double OverallMargin;
    double ISGSDiffMeanAvg;
    double ISGSDiffMeanMax;
    float IrisPupilGSDiff;
    float PupilIrisDiameterRatio;

    int pointresponse;
    int maxpointresponse;
};

/**
 * Returns the overall quality of the image.
 *
 * @return an integer indicating the overall quality of the image.
 */
/* inline */ int MFilter::getOverallQuality() { return OverallQuality; }

/**
 *
 * @return
 */
/* inline */ int MFilter::getDefocusScore() { return DefocusScore; }

/**
 *
 * @return
 */
/* inline */ int MFilter::getContrastScore() { return ContrastScore; }

/**
 * Returns the radius of the iris in the image.
 *
 * @return the radius in pixels.
 */
/* inline */ int MFilter::getIrisRadius() { return IrisRadius; }

/**
 * Returns the x coordinate of the center of the iris.
 *
 * @return the x coordinate in pixels
 */
/* inline */ int MFilter::getIrisCenterX() { return IrisCenterX; }

/**
 * Returns the y coordinate of the center of the iris.
 *
 * @return the y coordinate in pixels.
 */
/* inline */ int MFilter::getIrisCenterY() { return IrisCenterY; }

/**
 *
 * @return
 */
/* inline */ int MFilter::getUsableIrisAreaPercent()
{
    return UsableIrisAreaPercent;
}

/**
 *
 * @return
 */
/* inline */ float MFilter::getNDefocus() { return NDefocus; }

/**
 *
 * @return
 */
/* inline */ float MFilter::getNContrast() { return NContrast; }

/**
 *
 * @return
 */
/* inline */ float MFilter::getNISGSMean() { return NISGSMean; }

/**
 *
 * @return
 */
/* inline */ float MFilter::getIrisPupilGSDiff() { return IrisPupilGSDiff; }

/**
 *
 * @return
 */
/* inline */ float MFilter::getNIPGSDiff() { return NIPGSDiff; }

/**
 *
 * @return
 */
/* inline */ float MFilter::getNIrisID() { return NIrisID; }

/**
 *
 * @return
 */
/* inline */ float MFilter::getNIrisVis() { return NIrisVis; }

/**
 *
 * @return
 */
/* inline */ double MFilter::getISGSDiffMeanAvg() { return ISGSDiffMeanAvg; }

#endif
