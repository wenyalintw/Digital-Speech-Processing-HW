#include <iostream>
#include "hmm.h"
#include <fstream>

using namespace std;

template <class T>
T **make2Darray(int nr,int nc)
{ int size,i; T **m;
    m=(T **)calloc(nr, sizeof(T *));
    for(i=0;i<nr;i++){m[i]=(T *)calloc(nc, sizeof(T));}
    return m;
}

template <class T>
void free2Darray(T **m, int nr,int nc)
{   int i;
    for(i=0;i<nr;i++){free(m[i]);}
    free(m);
}

template <class T>
T ***make3Darray(int nr,int nc, int z)
{   T ***m;
    m=(T ***)calloc(nr, sizeof(T **));
    for(int i=0;i<nr;i++){m[i]=(T **)calloc(nc, sizeof(T *));}
    for(int i=0;i<nr;i++){
        for(int j=0;j<nc;j++){
            m[i][j]=(T *)calloc(z, sizeof(T));
        }
    }
    return m;
}

template <class T>
void free3Darray(T*** m, int nr,int nc, int z)
{
    for(int i=0;i<nr;i++){
        for(int j=0;j<nc;j++){
            free(m[i][j]);
        }
    }
    for(int i=0;i<nr;i++){free(m[i]);}
    free(m);
}


int getseqlen(char* trainingsequence){
    FILE *fp = open_or_die( trainingsequence, "r");
    char token[MAX_LINE] = "";
    // fscanf讀取fp中的string（%s）直到空白換行tab等截止，存到token中，回傳成功配對數目，所以只用一個%s要嘛1要嘛0
    fscanf( fp, "%s", token );
    int i=0;
    for (i=0; token[i]!='\0'; i++){}
    fclose(fp);
    return i;
}


struct load_result{
    int samples_count;
    int** train_data;
};

load_result loadSeq(char* trainingsequence, int T){
    int cols = T;

    fstream fin;
    fin.open(trainingsequence,ios::in);
    fin.seekg(0,ios::end);
    int char_count = fin.tellg();
    fin.seekg(0, ios::beg);

    int samples_count = char_count/(cols+1);
    int** train_data = make2Darray<int>(samples_count,cols);


    char linebuffer[100]="";                 //讀檔時暫存用
    int i=0;
    while(fin.getline(linebuffer, sizeof(linebuffer), '\n')) {
        for (int j = 0; j < 50; j++) {
            train_data[i][j] = linebuffer[j]-'A';
        }
        i++;
    }
    // 現在在train_data中，A~F分別對應數字0~5
    load_result result = {samples_count, train_data};
    return result;
}


void calcAlpha(double** alpha, double* pi, double** a, double** b, int N, int T, int** train_data, int current_seq){
    // Initial alpha
    // 初始化alpha first col (alpha: N * T ，N代表state，T為Discrete time)
    for (int i=0; i<N; i++){
        alpha[i][0] = pi[i] *  b[ train_data[current_seq][0] ][i];
    }

    // Induction alpha
    // 對所有t時刻的state i，下一時刻到j的機率總和，就是下一時刻的pi的概念，這時再去乘上observ噴出的值屬於state j的機率，即alpha[j][t+1]
    // 一路induct到alpha滿為止，所以要induct (N)*(T-1)次，每一次中所有i*a_ij要總合，所以會有三層loop
    for (int t=1; t<T; t++){
        for (int j=0; j<N; j++){
            // alpha[][]原始元素都是0了，直接先把i*a_ij累加上去，最後再乘上observ即可
            for (int i=0; i<N; i++){
                alpha[j][t] += alpha[i][t-1] * a[i][j];
            }
            alpha[j][t] *= b[ train_data[current_seq][t] ][j];
        }
    }
}

void calcBeta(double** beta, double* pi, double** a, double** b, int N, int T, int** train_data, int current_seq){
    // Initial alpha
    // 最後一cols全部都=1
    for (int i=0; i<N; i++){
        beta[i][T-1] = 1.;
    }
    // Induction
    // 3-loops，beta[i][t] = sum(所有state) (a_ij * observ@t+1 * beta[j][t+1 ])
    // 意思也就是: 從現在時間t我在state i，到t+1的時候可能會到各種state j，把我在t+1時可以到的全部機率乘起來即是
    // 每個state都要算，所以是三層
    for (int t=T-2; t>=0; t--){
        for (int i=0; i<N; i++){
            for (int j=0; j<N; j++){
                beta[i][t] += a[i][j] * b[ train_data[current_seq][t+1]][j] * beta[j][t+1];
            }
        }
    }
}

void calcGamma(double** gamma, double** alpha, double** beta, int N, int T){
    // 分母的sigma要一層loop，對每個t都做一次，共兩層
    for (int t=0; t<T; t++){
        double sum=0.;
        for(int i=0; i<N; i++){
            sum += (alpha[i][t] * beta[i][t]);
        }
        for(int i=0; i<N; i++){
            gamma[i][t] = (alpha[i][t] * beta[i][t] / sum);
        }
    }
}

void calcEpsilon(double*** epsilon, double** alpha, double** beta, double** a, double** b, int** train_data, int N, int T, int current_seq){

    for(int t=0; t<T-1; t++){
        double sum=0.;
        for(int i=0; i<N; i++){
            for(int j=0; j<N; j++){
                sum += alpha[i][t] * a[i][j] * b[ train_data[current_seq][t + 1] ][j] *  beta[j][t+1];
            }
        }
        for(int i=0; i<N; i++){
            for(int j=0; j<N; j++){
                epsilon[i][j][t] = (alpha[i][t] * a[i][j] * b[ train_data[current_seq][t + 1] ][j] *  beta[j][t+1] / sum );
            }
        }
    }
}

void accumulate(double** gamma, double*** epsilon, double* pi_num, double** a_num, double* a_den,  double** b_num,  double* b_den, int** train_data,int N, int T, int current_seq){

    for (int i = 0; i < N; i++)
        pi_num[i] += gamma[i][0];


    for (int i = 0; i < N; i++)
        for(int t = 0; t < T - 1; t++)
            a_den[i] += gamma[i][t];
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            for(int t = 0; t < T - 1; t++)
                a_num[i][j] += epsilon[i][j][t];


    for (int j = 0; j < N; j++)
        for(int t = 0; t < T; t++)
            b_den[j] += gamma[j][t];
    for (int j = 0; j < N; j++)
        for (int t = 0; t < T; t++)
            b_num[ train_data[current_seq][t] ][j] += gamma[j][t];
}

void updateParas(double* pi, double** a, double** b, double* pi_num, double** a_num, double* a_den,  double** b_num,  double* b_den, int N, int samples_count, int M){

    // update pi
    for (int i = 0; i<N; i++){
        pi[i] = pi_num[i] / samples_count;
    }

    // update a
    for (int i=0; i<N; i++){
        for (int j=0; j<N; j++){
            a[i][j] = a_num[i][j] / a_den[i];
        }
    }

    // update b
    for (int i=0; i<M; i++){
        for (int j=0; j<N; j++){
            b[i][j] = b_num[i][j] / b_den[j];
        }
    }
}


int main(int argc, char **argv) {


    setbuf(stdout, 0);
    /* 讀入command line argument，檢測無問題就存入變數 */
    if (argc!=5){
        cerr << "Usage: " << argv[0] << "<iterations> <model_init.txt> <seq_model_0X.txt> <model_0X.txt>";
        exit (EXIT_FAILURE);
    }
    long int iterations = strtol(argv[1], nullptr, 10);
    char initialmodel[40];// = strdup(argv[2]);
    char trainingsequence[40];// = strdup(argv[3]);
    char outfile[40];// = strdup(argv[4]);
    strcpy(initialmodel, argv[2]);
    strcpy(trainingsequence, argv[3]);
    strcpy(outfile, argv[4]);


    // 建立一個HMM model
    // struct沒有constructor，所以這邊用pointer是自找麻煩，後續pass進去的時候還沒allocate memory，直接create object實在
    HMM hmm;
    // 拿initial model去初始化它
    loadHMM(&hmm, initialmodel);

    //const int T = 50;
    const int T = getseqlen(trainingsequence);
    const int N = hmm.state_num;
    const int M = hmm.observ_num;

    load_result temp = loadSeq(trainingsequence, T);
    int** train_data = temp.train_data;
    int samples_count = temp.samples_count;
    double* pi = hmm.initial;

    // a,b
    double** a = make2Darray<double>(MAX_STATE, MAX_STATE);
    double** b = make2Darray<double>(MAX_OBSERV, MAX_STATE);
    double (*a_temp)[MAX_STATE] = hmm.transition;
    double (*b_temp)[MAX_STATE] = hmm.observation;
    for(int i=0; i<MAX_STATE; i++){
        for(int j=0; j<MAX_STATE; j++){a[i][j] = a_temp[i][j];}
    }
    for(int i=0; i<MAX_OBSERV; i++){
        for(int j=0; j<MAX_STATE; j++){b[i][j] = b_temp[i][j];}
    }

    // 開始計算
    for (int iter = 0; iter < iterations; iter++) {
        // 單次iter內要做的事
        auto pi_num = (double *) calloc(N, sizeof(double));
        double **a_num = make2Darray<double>(N, N);
        auto a_den = (double *) calloc(N, sizeof(double));
        double **b_num = make2Darray<double>(M, N);
        auto b_den = (double *) calloc(N, sizeof(double));

        for (int current_seq = 0; current_seq < samples_count; current_seq++) {
            double** alpha = make2Darray<double>(N,T);
            double** beta = make2Darray<double>(N,T);
            double** gamma = make2Darray<double>(N,T);
            double*** epsilon = make3Darray<double>(N,N,T-1);
            calcAlpha(alpha, pi, a, b, N, T, train_data, current_seq);
            calcBeta(beta, pi, a, b, N, T, train_data, current_seq);
            calcEpsilon(epsilon, alpha, beta, a, b, train_data, N, T, current_seq);
            calcGamma(gamma, alpha, beta, N, T);
            accumulate(gamma, epsilon, pi_num, a_num, a_den, b_num, b_den, train_data, N, T, current_seq);

            free2Darray<double>(alpha, N, T);free2Darray<double>(beta, N, T);free2Darray<double>(gamma, N, T);free3Darray<double>(epsilon, N, N, T-1);
            //free(alpha); free(beta); free(gamma); free(epsilon);
        }
        // 每次iter都要update para
        updateParas(pi, a, b, pi_num, a_num, a_den, b_num, b_den, N, samples_count, M);
        // update用參數重置
        free(pi_num); free(a_den);free(b_den);
        free2Darray<double>(a_num, N, N); free2Darray<double>(b_num, M, N);

    }

    // 把a,b設回hmm
    for(int i=0; i<MAX_STATE; i++){
        for(int j=0; j<MAX_STATE; j++){a_temp[i][j] = a[i][j];}
    }
    for(int i=0; i<MAX_OBSERV; i++){
        for(int j=0; j<MAX_STATE; j++){b_temp[i][j] = b[i][j];}
    }

    // 輸出output檔
    FILE *fout = fopen(argv[4], "wb");
    dumpHMM(fout, &hmm);
    fclose(fout);

    //把為了argv，alloc的變數記憶體釋放
    //竟然沒free zzz
    free(hmm.model_name);
    free2Darray<int>(train_data, samples_count, T); free2Darray<double>(a, MAX_STATE, MAX_STATE); free2Darray<double>(b, MAX_OBSERV, MAX_STATE);
    return 0;
}