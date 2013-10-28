#include <unistd.h>
#include <stdlib.h>
#include <opencv/cv.h>
#include <opencv/highgui.h>
#include <fstream>
#include <getopt.h>
#include <assert.h>
#include <string>
#include <sys/time.h>
#include <time.h>

using namespace cv;
using namespace std;

#include <cilk/cilk.h>
#include <cilk/cilk_api.h>

extern "C" {
#include "gaussian_conv_vyv.h"
}

void adapthisteq(const Mat &in, Mat &out, float regularizer);
void local_statistics(Mat &image_in, int windowsize, Mat &mean, Mat &var, Mat deciles[10]);

void dog_2_50(const Mat &in, Mat &out)
{
    Mat blur_2, blur_50, diff;
    cilk_spawn GaussianBlur(in, blur_2, Size(0, 0), 2);
    cilk_spawn GaussianBlur(in, blur_50, Size(0, 0), 50);
    cilk_sync;
    subtract(blur_2, blur_50, diff, noArray(), CV_16S);
    diff.convertTo(out, CV_8U, 0.5, 128);
}

void vyv_blur(const Mat &in, Mat &out, double sigma)
{
    vyv_coeffs c;
    int K = 3;
    double tol = 1e-2;

    out.create(in.size(), in.type());
    vyv_precomp(&c, sigma, K, tol);

    /* Blur rows, store in out */
    cilk_for (int i = 0; i < in.rows; i++) {
        Mat tmp;
        in.row(i).convertTo(tmp, CV_32F);
        float *p = (float *) tmp.ptr();
        vyv_gaussian_conv(c, p, p, tmp.cols, 1);
        tmp.convertTo(out.row(i), in.type());
    }

    /* Blur cols in place */
    cilk_for (int i = 0; i < out.cols; i++) {
        Mat tmp;
        out.col(i).convertTo(tmp, CV_32F);
        float *p = (float *) tmp.ptr();
        vyv_gaussian_conv(c, p, p, tmp.rows, 1);
        tmp.convertTo(out.col(i), out.type());
    }
}
    
    

void dog_2_50_vyv(const Mat &in, Mat &out)
{
    out.create(in.size(), in.type());
 
    Mat blur_2, blur_50, diff;
    vyv_blur(in, blur_2, 2.0);
    vyv_blur(in, blur_50, 50.0);
    subtract(blur_2, blur_50, diff, noArray(), CV_16S);
    diff.convertTo(out, CV_8U, 0.5, 128);
}

int main(int argc, char **argv)
{
    int iters = 200;

    /* Read input, convert to grayscale */
    Mat im;
    im = imread(argv[1], 0);
    im.convertTo(im, CV_8U);


    clock_t cpu_begin, cpu_end;
    struct timeval wall_start, wall_end;
    double cpu_secs, wall_secs;
    
    
    cpu_begin = clock();
    gettimeofday(&wall_start, NULL);
    for (int i = 0; i < iters; i++) {
        Mat tmp;
        dog_2_50(im, tmp);
    }
    cpu_end = clock();
    gettimeofday(&wall_end, NULL);
    cpu_secs = (double)(cpu_end - cpu_begin) / CLOCKS_PER_SEC;
    wall_secs = wall_end.tv_sec + 1e-6 * wall_end.tv_usec - \
        (wall_start.tv_sec + 1e-6 * wall_start.tv_usec);

    cout << iters << " iters cv2:" << endl;
    cout << "  wall: " << wall_secs << endl;
    cout << "  CPU secs: " << cpu_secs << endl;
    cout << "  CPU/wall: " << cpu_secs / wall_secs << endl;

    cpu_begin = clock();
    gettimeofday(&wall_start, NULL);
    for (int i = 0; i < iters; i++) {
        Mat tmp;
        dog_2_50_vyv(im, tmp);
    }
    cpu_end = clock();
    gettimeofday(&wall_end, NULL);
    cpu_secs = (double)(cpu_end - cpu_begin) / CLOCKS_PER_SEC;
    wall_secs = wall_end.tv_sec + 1e-6 * wall_end.tv_usec - \
        (wall_start.tv_sec + 1e-6 * wall_start.tv_usec);
    cout << iters << " iters vyv:" << endl;
    cout << "  wall: " << wall_secs << endl;
    cout << "  CPU secs: " << cpu_secs << endl;
    cout << "  CPU/wall: " << cpu_secs / wall_secs << endl;


}
