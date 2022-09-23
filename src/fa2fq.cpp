#include "HaploCart.h"
#include "miscfunc.h"
#include<random>


void Haplocart::write_fq_read(auto & dummyFASTQFile, int offset, const int window_size, const string &fastaseq, char dummyqualscore) {
        string seq_to_write, qual_to_write;
        for (const char base : fastaseq.substr(min(offset, int(fastaseq.size())), window_size)) {
            if (base != 'N') {
                seq_to_write += base;
                qual_to_write += dummyqualscore;
                             }
            else {
                  seq_to_write += 'A';
                  qual_to_write += '!';
                 }
                                                               }

        dummyFASTQFile << '@' << random_string(7) << '\n';
        dummyFASTQFile << seq_to_write;
        dummyFASTQFile << "\n+\n";
        dummyFASTQFile << qual_to_write;
        dummyFASTQFile << '\n';
        offset += window_size;
                     }


const string Haplocart::fa2fq(const string & fastaseq, const char & dummyqualscore, const string & tmpdir) {

    // Given a consensus sequence and a dummy quality score, convert to FASTQ

    // https://stackoverflow.com/questions/7560114/random-number-c-in-some-range

    const string prefix = tmpdir;
    const string tempfqfilename= prefix+random_string(7);
    ofstream dummyFASTQFile;
    dummyFASTQFile.open(tempfqfilename);
    int window_size = ceil(fastaseq.size()/100);
    int offset = 0;
    string seq_to_write, qual_to_write;

     for (int i = 0; i < 101; ++i) {
        Haplocart::write_fq_read(dummyFASTQFile, offset, window_size, fastaseq, dummyqualscore);
        offset += 100;
                                   }

     for (int i = 1; i < 101; ++i) {
        Haplocart::write_fq_read(dummyFASTQFile, offset, window_size, fastaseq, dummyqualscore);
        offset += 100;
                                   }

     dummyFASTQFile.close();
     return tempfqfilename;
                                                                                     }


