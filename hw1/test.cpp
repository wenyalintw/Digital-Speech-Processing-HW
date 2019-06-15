//
// Created by 温雅 on 2019-03-28.
//

// https://zh.wikipedia.org/wiki/%E7%BB%B4%E7%89%B9%E6%AF%94%E7%AE%97%E6%B3%95
// 尋找最有可能產生觀測事件序列的維特比路徑，每個model都得到該個路徑的機率，比較每個model的這個機率
#include "hmm.h"
#include <iostream>
#include <fstream>
#include <iomanip>

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


struct load_result{
    int samples_count;
    int** test_data;
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


struct viterbi_result{
    char* model_name;
    double value;
};


viterbi_result viterbi(HMM* hmmlist, int** test_data, int N, int T, int current_seq){
    double P_final = 0.;
    int c_final;
    for (int c=0; c<5; c++) {
        double** delta =  make2Darray<double>(N,T);

        // Initialization
        for (int i = 0; i < N; i++) {
            delta[i][0] = hmmlist[c].initial[i] * hmmlist[c].observation[test_data[current_seq][0]][i];
        }
        // Recursion
        for (int t = 1; t < T; t++) {
            for (int j = 0; j < N; j++) {
                double max = 0.;
                for (int i = 0; i < N; i++) {
                    double temp = delta[i][t - 1] * hmmlist[c].transition[i][j];
                    max = temp > max? temp:max;
                }
                delta[j][t] = max * hmmlist[c].observation[test_data[current_seq][t]][j];
            }
        }
        // Termination
        double P = 0.;
        for (int i = 0; i < N; i++) {
            P = delta[i][T-1] > P ? delta[i][T-1] : P;
        }

        if (P>P_final){
            P_final = P;
            c_final = c;
        }

        free2Darray(delta, N, T);
    }

    viterbi_result result = {hmmlist[c_final].model_name, P_final};
    return result;
}

int main(int argc, char **argv) {


    setbuf(stdout, 0);
    /* 讀入command line argument，檢測無問題就存入變數 */
    if (argc!=4){
        cerr << "Usage: " << argv[0] << "<modellist.txt> <testing_data.txt> <result.txt>";
        exit (EXIT_FAILURE);
    }
    char modellist[40];
    char testingsequence[40];
    char outfile[40];
    strcpy(modellist, argv[1]);
    strcpy(testingsequence, argv[2]);
    strcpy(outfile, argv[3]);


    // load 5個model
    HMM hmmlist[5];
    load_models( modellist, hmmlist, 5 );

    // 讀testing data
    // 先確定單行長度
    const int T = getseqlen(testingsequence);
    // load進來
    load_result temp = loadSeq(testingsequence, T);
    int** test_data = temp.test_data;
    int samples_count = temp.samples_count;
    cout << endl;

    // N
    const int N = hmmlist[0].state_num;

    fstream fout;
    fout.open(outfile, ios::out | ios::trunc);
    for (int current_seq = 0; current_seq < samples_count; current_seq++) {
        viterbi_result seq_result = viterbi(hmmlist, test_data, N, T, current_seq);
        fout << seq_result.model_name << " " << scientific << seq_result.value << endl ;
    }
    fout.close();

    free2Darray(test_data, samples_count, T);
    for(int i=0; i<5; i++){
        free(hmmlist[i].model_name);
    }

    return 0;
}