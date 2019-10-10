/* ep2d_simple_15.cpp */
#include "mex.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <cmath>
#include <thread>
#include <vector>
#include <valarray>

static const double CONST_ITR_SCALE = 0.50;
static const int CONST_FILTER_SIZE = 15;
static const int CONST_FILTER_SIZE_HALF = 7;
using namespace std;

static inline
void filter_gray(const int j, const int M, const double *temporary_img_yoko, const double *temporary_img_tate, const double threshold, double *output,const double *input_img);
static inline
void filter_color(const int j, const int M, const int N, const double *temporary_img_yoko, const double *temporary_img_tate, const double threshold, double *output,const double *input_img);


void ep2d_simple_gray(double *output, double *input, const int M, const int N, double threshold, const int iter_max){
    
    int i,j,k,l;
    int index = 0;
    int iter;
    double *temporary_img_yoko,*temporary_img_tate,*input_img;
    double computed_sum=0.0,filter_sum=0.0;
    const int TempM = M+(CONST_FILTER_SIZE-1);
    
    temporary_img_yoko = (double *)mxCalloc((M+(CONST_FILTER_SIZE-1))*(N+(CONST_FILTER_SIZE-1)),sizeof(double));
    temporary_img_tate = (double *)mxCalloc((M+(CONST_FILTER_SIZE-1))*(N+(CONST_FILTER_SIZE-1)),sizeof(double));
    input_img = (double *)mxCalloc((M+(CONST_FILTER_SIZE-1))*(N+(CONST_FILTER_SIZE-1)),sizeof(double));
    memcpy(output, input, sizeof(double) * M*N);
    
    for (iter=0;iter<iter_max;++iter){
        
        // initialization
        for (j=CONST_FILTER_SIZE_HALF; j<N+CONST_FILTER_SIZE_HALF; ++j){
            memcpy(&input_img[CONST_FILTER_SIZE_HALF + j*TempM],&output[(j-CONST_FILTER_SIZE_HALF)*M],sizeof(double) * M);
        }
        
        // temporaryの作成
        for (j=1; j<N+CONST_FILTER_SIZE_HALF; ++j){
            for (i=1; i<M+CONST_FILTER_SIZE_HALF; ++i){
                index = i+j*TempM;
                temporary_img_yoko[index] = fabs(input_img[i+j*TempM]-input_img[index-TempM]) + temporary_img_yoko[index-TempM];
                temporary_img_tate[index] = fabs(input_img[i+j*TempM]-input_img[index-1]) + temporary_img_tate[index-1];
            }
        }
        
        // filtering
        int MaxThread = max(thread::hardware_concurrency(), 1u);
        vector<thread> worker;
        for (int ThreadNum = 0; ThreadNum < MaxThread; ThreadNum++) {
            worker.emplace_back([&](int id) {
                int r0 = N/MaxThread * id + min(N%MaxThread, id);
                int r1 = N/MaxThread * (id+1) + min(N%MaxThread, id+1);
                for (int j = r0; j < r1; ++j) {   
                    filter_gray(j,M,temporary_img_yoko,temporary_img_tate,threshold,output,input_img);
                }
            }, ThreadNum);}
        for (auto& t : worker) t.join();
        
        if (iter<iter_max){
            threshold = threshold*CONST_ITR_SCALE;
        }
        
    }
    
    
    // end process
    mxFree(temporary_img_yoko);
    mxFree(temporary_img_tate);
    mxFree(input_img);
    
}

static inline
void filter_gray(const int j, const int M, const double *temporary_img_yoko, const double *temporary_img_tate, const double threshold, double *output,const double *input_img){
    double filter_rcon_yoko[CONST_FILTER_SIZE*CONST_FILTER_SIZE],filter_rcon_tate[CONST_FILTER_SIZE*CONST_FILTER_SIZE];
    double computed_sum,filter_sum;
    int index;
    const int TempM = M+(CONST_FILTER_SIZE-1);

    
    
    for (int i=0; i<M; ++i){
        // reconstructed kernel
        for (int l=0;l<CONST_FILTER_SIZE;++l){
            for (int k=0;k<CONST_FILTER_SIZE;++k){
                index = (i+k)+(j+l)*TempM;
                filter_rcon_yoko[k+l*CONST_FILTER_SIZE] = fabs(temporary_img_yoko[index]-temporary_img_yoko[(i+k)+(j+CONST_FILTER_SIZE_HALF)*TempM]);
                filter_rcon_tate[k+l*CONST_FILTER_SIZE] = fabs(temporary_img_tate[index]-temporary_img_tate[index-k+CONST_FILTER_SIZE_HALF]);
            }
        }
        
        // compute&apply filter
        computed_sum = 0.0;
        filter_sum = 0.0;
        for (int l=0;l<CONST_FILTER_SIZE;++l){
            for (int k=0;k<CONST_FILTER_SIZE;++k){
                index = k+l*CONST_FILTER_SIZE;
                if (fmin(filter_rcon_yoko[index] + filter_rcon_tate[k+(CONST_FILTER_SIZE_HALF*CONST_FILTER_SIZE)],filter_rcon_tate[index] + filter_rcon_yoko[index-k+CONST_FILTER_SIZE_HALF]) < threshold){
                    computed_sum += input_img[(i+k)+(j+l)*TempM];
                    filter_sum++;
                }
            }
        }
        output[i+j*M] = computed_sum/filter_sum;
    }
    

}




void ep2d_simple_color(double *output, double *input, const int M, const int N, double threshold, const int iter_max){
    
    int index = 0;
    const int TempLayerSize = (M+(CONST_FILTER_SIZE-1))*(N+(CONST_FILTER_SIZE-1));
    const int TempM = (M+(CONST_FILTER_SIZE-1));
    int iter = 0;
    double *temporary_img_yoko,*temporary_img_tate,*input_img;
    double *filter_rcon_yoko,*filter_rcon_tate;
    
    temporary_img_yoko = (double *)calloc(TempLayerSize,sizeof(double));
    temporary_img_tate = (double *)calloc(TempLayerSize,sizeof(double));
    input_img = (double *)calloc(TempLayerSize*3,sizeof(double));
    memcpy(output, input, sizeof(double) * M*N*3);
    
    for (iter=0;iter<iter_max;iter++){
        
        // initialization
        for (int j=CONST_FILTER_SIZE_HALF; j<N+CONST_FILTER_SIZE_HALF; j=j+1){
            memcpy(&input_img[CONST_FILTER_SIZE_HALF + j*TempM],&output[(j-CONST_FILTER_SIZE_HALF)*M],sizeof(double) * M);
            memcpy(&input_img[CONST_FILTER_SIZE_HALF + j*TempM+TempLayerSize],&output[(j-CONST_FILTER_SIZE_HALF+N)*M],sizeof(double) * M);
            memcpy(&input_img[CONST_FILTER_SIZE_HALF + j*TempM+2*TempLayerSize],&output[(j-CONST_FILTER_SIZE_HALF+2*N)*M],sizeof(double) * M);
        }
        
        // temporaryの作成
        index = 0;
        for (int j=1; j<N+CONST_FILTER_SIZE_HALF; ++j){
            index = j*TempM;
            for (int i=1; i<M+CONST_FILTER_SIZE_HALF; ++i){
                ++index;
                temporary_img_yoko[index] = fabs(input_img[index]-input_img[index-TempM])
                + fabs(input_img[index+TempLayerSize]-input_img[index-TempM+TempLayerSize])
                + fabs(input_img[index+TempLayerSize+TempLayerSize]-input_img[index-TempM+TempLayerSize+TempLayerSize])
                + temporary_img_yoko[index-TempM];
                temporary_img_tate[index] = fabs(input_img[index]-input_img[index-1])
                + fabs(input_img[index+TempLayerSize]-input_img[index-1+TempLayerSize])
                + fabs(input_img[index+TempLayerSize+TempLayerSize]-input_img[index-1+TempLayerSize+TempLayerSize])
                + temporary_img_tate[index-1];
            }
        }
        
        // filtering
        int MaxThread = max(thread::hardware_concurrency(), 1u);
        vector<thread> worker;
        for (int ThreadNum = 0; ThreadNum < MaxThread; ThreadNum++) {
            worker.emplace_back([&](int id) {
                int r0 = N/MaxThread * id + min(N%MaxThread, id);
                int r1 = N/MaxThread * (id+1) + min(N%MaxThread, id+1);
                for (int j = r0; j < r1; ++j) {   
                    filter_color(j,M,N,temporary_img_yoko,temporary_img_tate,threshold,output,input_img);
                }
            }, ThreadNum);}
        for (auto& t : worker) t.join();
        
        
        if (iter<iter_max){
            threshold = threshold*CONST_ITR_SCALE;
        }
        
    }
    
    
    
    // end process
    free(temporary_img_yoko);
    free(temporary_img_tate);
    free(input_img);
    
}


static inline
void filter_color(const int j, const int M, const int N, const double *temporary_img_yoko, const double *temporary_img_tate, const double threshold, double *output,const double *input_img){
    double filter_rcon_yoko[CONST_FILTER_SIZE*CONST_FILTER_SIZE],filter_rcon_tate[CONST_FILTER_SIZE*CONST_FILTER_SIZE];
    valarray<double> filter_sumsum(0.0, 3);
    int index = 0;
    int index2 = 0;
    int id_count = 0;
    const int TempLayerSize = (M+(CONST_FILTER_SIZE-1))*(N+(CONST_FILTER_SIZE-1));
    const int TempM = (M+(CONST_FILTER_SIZE-1));
    
    
    for (int i=0; i<M; ++i){
                
                // reconstructed kernel
                for (int l=0;l<CONST_FILTER_SIZE;++l){
                    index = i+(j+l)*TempM;
                    for (int k=0;k<CONST_FILTER_SIZE;++k){
                        filter_rcon_yoko[k+l*CONST_FILTER_SIZE] = fabs(temporary_img_yoko[index]-temporary_img_yoko[(i+k)+(j+CONST_FILTER_SIZE_HALF)*TempM]);
                        filter_rcon_tate[k+l*CONST_FILTER_SIZE] = fabs(temporary_img_tate[index]-temporary_img_tate[index-k+CONST_FILTER_SIZE_HALF]);
                        ++index;
                    }
                }
                
                //compute&apply filter
                filter_sumsum = 0.0;
                index = 0;
                id_count = 0;
                for (int l=0;l<CONST_FILTER_SIZE;++l){
                    for (int k=0;k<CONST_FILTER_SIZE;++k){
                        if (fmin(filter_rcon_tate[index] + filter_rcon_yoko[index-k+CONST_FILTER_SIZE_HALF],filter_rcon_yoko[index] + filter_rcon_tate[k+(CONST_FILTER_SIZE_HALF*CONST_FILTER_SIZE)]) < threshold){
                            index2 = (i+k)+(j+l)*TempM;
                            ++id_count;
                            filter_sumsum[0] += input_img[index2];
                            filter_sumsum[1] += input_img[index2 + TempLayerSize];
                            filter_sumsum[2] += input_img[index2 + TempLayerSize + TempLayerSize];
                        }
                        ++index;
                    }
                }

                filter_sumsum /= id_count;
                index = i+j*M;
                output[index] = filter_sumsum[0];
                output[index + M*N] = filter_sumsum[1];
                output[index + 2*M*N] = filter_sumsum[2];
            }
    
}



void mexFunction( int Nreturned, mxArray *returned[], int Noperand, const mxArray *operand[] ){
    double *input, *output;
    const mwSize *input_size_temp;
    mwSize input_size[3];
    mwSize dim_count;
    
    
    /* 画像の次元数の判定 */
    dim_count = mxGetNumberOfDimensions(operand[0]);
    input_size_temp = mxGetDimensions(operand[0]);
    if (dim_count < 3){
        input_size[0] = input_size_temp[0];
        input_size[1] = input_size_temp[1];
        input_size[2] = 1;
    }else{
        input_size[0] = input_size_temp[0];
        input_size[1] = input_size_temp[1];
        input_size[2] = input_size_temp[2];
    }
    
    // 返り値用の領域を確保
    returned[0] = mxCreateNumericArray(dim_count,input_size, mxDOUBLE_CLASS, mxREAL);
    
    // Matlab側の変数のアドレスをC側の変数にコピーする
    input = mxGetPr(operand[0]);
    output = mxGetPr(returned[0]);
    
    //multithread
    //printf("Processor:%d\n",std::thread::hardware_concurrency());
    
    // 処理を実行
    if (input_size[2] == 1){
        ep2d_simple_gray(output,input,input_size[0],input_size[1],(double) *mxGetPr(operand[1]),(int) *mxGetPr(operand[2]));
    }else{
        ep2d_simple_color(output,input,input_size[0],input_size[1],(double) *mxGetPr(operand[1]),(int) *mxGetPr(operand[2]));
    }
}
