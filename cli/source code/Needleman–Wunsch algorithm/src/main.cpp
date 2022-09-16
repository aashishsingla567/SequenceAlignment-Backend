// CLI to execute needleman-wunsch algorithm
// on nucleotide sequences
// and show the Alignment

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "./includes/json.hpp" // for pasing and stringifying JSON


template <typename T>
using matrix = std::vector < std::vector <T> >;

// Scoring schema
struct Scoring {
    int match;
    int mismatch;
    int gap;
};

// default Scoring schema
inline Scoring default_scoring() {
    return {
         1, // match
        -1, // mismatch
         0  // gap
    };
}

struct Alignment {
    std::string seq1;
    std::string seq2;
    matrix <int> score;
public:
    int calc_score (Scoring scr) {
        int score = 0;
        for (int i = 0; i < seq1.size(); i++) {
            if (seq1[i] == seq2[i]) {
                score += scr.match;
            } else if (seq1[i] == '-' || seq2[i] == '-') {
                score += scr.gap;
            } else {
                score += scr.mismatch;
            }
        }
        return score;
    }
};

matrix <int> generateAlignmentMatrix (int64_t m, int64_t n, Scoring scr) {
    matrix <int> mat (m + 1, std::vector <int> (n + 1, 0));
    for (int64_t i = 1; i <= m; ++i) {
        mat[i][0] = mat[i - 1][0] + scr.gap;
    }
    for (int64_t j = 1; j <= n; ++j) {
        mat[0][j] = mat[0][j - 1] + scr.gap;
    }
    return mat;
}

Alignment alignment_algorithm (const std::string& seq1, const std::string& seq2, Scoring scr = default_scoring()) {
    matrix <int> mat = generateAlignmentMatrix (seq1.size(), seq2.size(), scr);

    // fill the matrix
    for (int64_t i = 1; i <= seq1.size(); ++i) {
        for (int64_t j = 1; j <= seq2.size(); ++j) {
            int match = mat[i - 1][j - 1] + (seq1[i - 1] == seq2[j - 1] ? scr.match : scr.mismatch);
            int del = mat[i - 1][j] + scr.gap;
            int ins = mat[i][j - 1] + scr.gap;
            mat[i][j] = std::max (match, std::max (del, ins));
        }
    }

    // tracebacks
    std::string al1, al2;
    // reserve string size
    al1.reserve (std::max (seq1.size(), seq2.size()));
    al2.reserve (std::max (seq1.size(), seq2.size()));

    int64_t i = seq1.size(), j = seq2.size();

    while (i > 0 && j > 0) {
        // traceback
        int score = mat[i][j];
        int diag = mat[i - 1][j - 1];
        int up = mat[i][j - 1];
        int left = mat[i - 1][j];

        if (seq1[i - 1] == seq2[j - 1] || diag >= up && diag >= left) {
            // move diagonally
            al1.push_back (seq1[i - 1]);
            al2.push_back (seq2[j - 1]);
            --i;
            --j;
        } else if (left >= up && left >= diag) {
            // move left
            al1.push_back (seq1[i - 1]);
            al2.push_back ('-');
            --i;
        } else if (up >= left && up >= diag) {
            // move up
            al1.push_back ('-');
            al2.push_back (seq2[j - 1]);
            --j;
        }

    }
    // copy remaining characters
    while (i > 0) {
        al1.push_back (seq1[i - 1]);
        al2.push_back ('-');
        --i;
    }
    while (j > 0) {
        al1.push_back ('-');
        al2.push_back (seq2[j - 1]);
        --j;
    }
    // reverse the strings
    std::reverse (al1.begin(), al1.end());
    std::reverse (al2.begin(), al2.end());
    return {al1, al2, mat};
}

void printResults (Alignment& al, int score) {

    // print results
    std::cout << "Results " << "\n\n";
    // print score
    std::cout << "Score: " << score << std::endl;

    // output the Alignment
    std::cout << "Seq1:: " << al.seq1 << std::endl;
    std::cout << "Seq2:: " << al.seq2 << std::endl;

    std::cout << "Matrix:: " << std::endl;
    // print matrix
    for (int i = 0; i < al.score.size(); i++) {
        for (int j = 0; j < al.score[i].size(); j++) {
            std::cout << al.score[i][j] << " \t";
        }
        std::cout << std::endl;
    }
}


int main (int argv, char* argc[]) {
    
    using namespace nlohmann;

    std::cout << "Reading files..." << std::endl;

    // take the input from a json file
    std::ifstream input (argc[1]);

    auto input_json = json::parse (input);
    
    std::string seq1 = input_json["seq1"];
    std::string seq2 = input_json["seq2"];
    Scoring scr_schema = {
        input_json["scoring_schema"]["match"],
        input_json["scoring_schema"]["mismatch"], 
        input_json["scoring_schema"]["gap"]
    };

    // print inputs
    std::cout << "Input:: " << std::endl;
    std::cout << input_json.dump (2) << '\n';
    std::cout << std::endl;
    // perform the alignment
    auto al = alignment_algorithm (seq1, seq2, scr_schema);

    int score = al.calc_score (scr_schema);

    // make a json from the results
    auto output_json = json::object({
        {"seq1", al.seq1},
        {"seq2", al.seq2},
        {"score", score},
        {"matrix", al.score}
    });
    std::cout << "Ouput:: " << std::endl;
    std::cout << output_json.dump (2) << std::endl;
    std::cout << std::endl;

    // write the json to a file
    std::ofstream output (argc[2]);
    output << output_json;

    return 0;
}