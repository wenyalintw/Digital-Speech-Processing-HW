//
// Created by �ථ on 2019-05-22.
//

#include <stdio.h>
#include "Ngram.h"
#include "VocabMap.h"

int main(int argc, char *argv[])
{
    setbuf(stdout, 0);
    Vocab vocab_bigram, vocab_zhuyin, vocab_big5;
    // Ngram�O�@��Language Model��class�A�Y�@��Vocab�Morder�]order�������ܹw�]�O3�^
    Ngram bigram_model(vocab_bigram, 2);

    // Ū�Jbigram.lm
    File bigram_lm("bigram.lm", "r" );
    bigram_model.read(bigram_lm);
    bigram_lm.close();

    // Ū�Jmap
    VocabMap zhuyin2big5_map(vocab_zhuyin, vocab_big5);
    File ZhuYin_Big5_map("ZhuYin-Big5.map", "r");
    zhuyin2big5_map.read(ZhuYin_Big5_map);
    ZhuYin_Big5_map.close();

    File testdata(argv[1], "r");
//    File testdata("testdata/1.txt", "r");
    // �i��B�z
    const int T = 200; //  �@��sentence�̦h�i�઺char��
    const int N = 1000; //  �@�Ӫ`���̦h�i������쪺Big5��
    char* current_line = nullptr;
    while((current_line = testdata.getline()) != nullptr){

        // �NŪ�J����eline�s��sentence���A�å[�W<s>�B<\s>
        VocabString sentence[T] = {};
        uint words_count = Vocab::parseWords(current_line, sentence+1, T); //+1��[0]����l�i�H��<s>
        sentence[0] = Vocab_SentStart;
        sentence[words_count+1] = Vocab_SentEnd;
        // ��scount�A�]���[�W�F�Y��
        words_count+=2;

        // �}�l��sentence�i��B�z�]Viterbi�^
        // LogP�N�Ofloat�A�]�����prob�Afunction�O�^��log���ȡA�o�˼g���n��{
        LogP probability[T][N] = {};
        // VocabIndex�N�Ounsigned int
        VocabIndex backtracking[T][N] = {};
        VocabIndex big5_index[T][N] = {};
        uint big5_count[T] = {};

        // Prob�N�Odouble
        Prob p_place_holder; VocabIndex current_index;
        // sentence[0]�û����O"<s>"�A��"<s>"��index�û����O1�A�w�q�_�l�I
        VocabMapIter iter(zhuyin2big5_map, 1);
        // next�|��scurrent_index�A�@�����ܷ�e��index
        iter.next(current_index, p_place_holder);
        probability[0][0] = -99; // -99��<s>�b.lm����LogP
        backtracking[0][0] = 0;
        big5_index[0][0] = current_index;
        big5_count[0]=1;

        VocabIndex vocab_unknown_index = vocab_bigram.getIndex(Vocab_Unknown);

        for (int t=1; t<words_count; t++){
            VocabMapIter zhuyin_row(zhuyin2big5_map, vocab_zhuyin.getIndex(sentence[t]));
            for (int j=0; zhuyin_row.next(current_index, p_place_holder); j++){
                big5_index[t][j] = current_index; big5_count[t]++;
                LogP bigram_max_prob = LogP_Zero;
                VocabIndex bigram_index1 = vocab_bigram.getIndex(vocab_big5.getWord(current_index));
                if (bigram_index1 == Vocab_None){ bigram_index1 = vocab_unknown_index; }
                for (int i=0; i<big5_count[t-1]; i++){
                    VocabIndex bigram_index2 = vocab_bigram.getIndex(vocab_big5.getWord(big5_index[t-1][i]));
                    if (bigram_index2 == Vocab_None){ bigram_index2 = vocab_unknown_index; }
                    VocabIndex context[] = {bigram_index2, Vocab_None};
                    LogP bigram_prob = bigram_model.wordProb(bigram_index1, context);
                    // �]�O�bLog scale�A���Ƭ��[
                    bigram_prob += probability[t-1][i];
                    if(bigram_prob > bigram_max_prob) {
                        bigram_max_prob = bigram_prob;
                        backtracking[t][j] = i;
                    }
                }
                probability[t][j] = bigram_max_prob;
            }
        }

        LogP path_max_prob = LogP_Zero;
        int final_index = 0;
        for (int j=0; j<big5_count[words_count-1]; j++){
            if (probability[words_count-1][j] > path_max_prob){
                path_max_prob = probability[words_count-1][j];
                final_index = j;
            }
        }

        int final_path_index[T];
        for (int t = words_count-1; t>=0; t--) {
            final_path_index[t] = final_index;
            final_index = backtracking[t][final_index];
        }

        for(int t = 0; t < words_count-1; t++){
            std::cout << vocab_big5.getWord(big5_index[t][ final_path_index[t] ] ) << " ";
        }
        std::cout << vocab_big5.getWord(big5_index[words_count-1][ final_path_index[words_count-1] ] ) << std::endl;
//        std::cout << "\b" << std::endl;
    }
    testdata.close();
    return 0;
}