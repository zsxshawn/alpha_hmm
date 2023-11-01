#include <array>
#include <cmath>
#include <fstream>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <map>
#include <tuple>
#include <queue>
#include <set>
#include <string>
#include <unordered_set>
#include <vector>
#include <unordered_map>
#include <filesystem>
#include <stack>



////Author: Aaron Yang, Sixiang
////Date: 2021-07-22 15:00:00
////Algorithm: alpha-HMM

using namespace std;
//ifstream fin ("test.in");
//ofstream fout ("test.out");
//// Emission matrix
//// Rows: State
//// Col: Observered
//double E[4][4] = {
//        {(double)1/4, (double)1/4, (double)1/4, (double)1/4},
//        {(double)1/4, (double)1/4, (double)1/4, (double)1/4},
//        {(double)1/4, (double)1/4, (double)1/4, (double)1/4},
//        {(double)1/4, (double)1/4, (double)1/4, (double)1/4},
//};
double E[4][4] = {
        {(double)0.22, (double)0.28, (double)0.29, (double)0.21},
        {(double)0.22, (double)0.28, (double)0.29, (double)0.21},
        {(double)0.22, (double)0.28, (double)0.29, (double)0.21},
        {(double)0.22, (double)0.28, (double)0.29, (double)0.21},
};
//// Transmission
//// Rows: from state
//// Cols: to state
const double Delta = 0.07;
const double Alpha = 0.2;
const double Beta = 0.015;
double T[4][4] = {
        {1 - Delta, Delta, 0, 0},
        {0, 0, 1, 0},
        {0, 0, 0, 1},
        {Alpha, Beta, 0, 1 - Alpha - Beta}
};


//// Influence probabilities based on observation
//// Rows: influencing observation
//// Cols: influenced observation
//// Order: ACGU
//double N[4][4] = {
//        {0, 0,  0,  0.72},
//        {0, 0,  0.96,  0},
//        {0, .96, 0, 0.32},
//        {0.72, 0, 0.32, 0},
//};
double N[4][4] = {
        {0, 0,  0,  0.125/0.22},
        {0, 0,  0.267/0.28,  0},
        {0,  0.416/0.29  ,0, 0.029/0.29},
        {0.142/0.21, 0, 0.021/0.21, 0},
};
//// Which states influence which states
//// Rows: influencing state
//// Cols: influenced state
bool allowedInfs[4][4]{
        {0, 1, 1, 1},
        {0, 0, 0, 0},
        {0, 0, 0, 0},
        {0, 0, 0, 0}
};

//// Probability of starting with each state
double startProb[4] = {1, 0, 0, 0};

//// Which end-states are allowed (1 is allowed, 0 is disallowed)
bool endState[4] = {1, 0, 0, 1};

//// Affiliations
//// Row # is affiliated with col #. Row# influence must be next to col#
bool aff[4][4] = {
        {0, 0, 0, 0},
        {0, 0, 0, 0},
        {0, 1, 0, 0},
        {0, 0, 1, 1},
};
//// Observation and state keys
vector<string> observationKey{"A", "C", "G", "U"};
vector<string> stateKey{"L", "P", "Q", "R"};

//// O is the input string
const int numStates = sizeof(T)/sizeof(T[0]);

std::vector<int> O;

//// Input string
void convertToNumbers(const std::string &input) {
    std::unordered_map<char, int> charMap = {
            {'A', 0},
            {'C', 1},
            {'G', 2},
            {'U', 3}
    };

    for (char c : input) {
        if (charMap.find(c) != charMap.end()) {
            O.push_back(charMap[c]);
        } else {
            std::cerr << "Invalid character " << c << " in the input string." << std::endl;
        }
    }
}
//// Check if influence is valid
//// influence: influence matrix
//// ingState: influencing state
//// edState: influenced state
//// return: true if valid, false if not
bool validInf(vector<vector<vector<int>>> &influence, int ingState, int edState, int ingStep, int edStep){
    //// Get current influence path
    //// Contains which steps are influencing or influenced

    //// Check if influenced is allowed. If not, return false.
    //// Use allowedInfs to check if influenced is allowed
    if (!allowedInfs[ingState][edState]){
        return false;
    }
    //// Check if step is already influencing some laterstep
    bool used = false;
    for (int i = 0; i < numStates; i++){
        for (int j = 0; j < numStates; j++){
            if (influence[i][j][ingStep] == 1){
                used = true;
                break;}
        }
    }
    return (!used);
}

//// Main function
std::string predictStructure(const std::string& sequence) {

    //// Input string
    //.........1.........2.........3.........4.........5.........6.........7.........8.........9.........0
    //"AGAAACUAGUUAAACUAAUAACACCGGAUUGUCAGACCGGAGUAACUGGUAAACAACCAGUGUUUCUUGCCA";

    //// Convert input string to numbers
    //cout << "Converting input string to numbers..." << endl;  // Debug statement
    convertToNumbers(sequence);
    //cout << "Input string: " << sequence << endl;  // Debug statement
    //// Number of steps, length of input string
    const int numSteps = O.size();
    //cout << "Number of steps: " << numSteps << endl;  // Debug statement
    //// Creates array to store which states are influenceable. If influence is possible, the influence will always be chosen.
    //// aka enforcing influences
    //// Rows: influenced state
    //// Cols: influencing state
    //// Equalvalent to influenceable = {false, true, true, true}
    bool influenceable[numStates];
    std::fill(influenceable, influenceable + numStates, false);
    //cout << "Influenceable states initialized." << endl;

    for (int i = 0; i < numStates; i++) {
        for (int j = 0; j < numStates; j++) {
            if (allowedInfs[i][j]) {
                influenceable[j] = true;
            }
        }
    }
    //cout << "Exited nested loop..." << endl;



    //// Probabilities
    ////m[current step][current state][influencing step]
    ////Size: (numSteps+1) * numStates * (numSteps+1)
    //// Initialize m function table
    //// Equivalent double m[numSteps+1][numStates][numSteps+1];
    vector<vector<vector<double>>> m(numSteps+1, vector<vector<double>>(numStates, vector<double>(numSteps+1, 0)));
   // cout << "m function table initialized." << endl;

    //// Influence matrix
    //// I[current step][current state][influencing step]
    //// Each entry of `I` has the following dimensions:
    //// 1. Current step
    //// 2. Current state
    //// 3. Influence step
    //// Each element of I[a][b][c][d] is a 3D vector:(vector<vector<vector<int>>>)

    vector<vector<vector<int>>> I[numSteps + 1][numStates][numSteps+1];

    //// Blank vector (value for each entry in `I`)
    //// Each entry of `blank` has the following dimensions:
    //// The value stored is a 3D vector, further broken down by:
    //// 1. Influencing state - which state is exerting the influence.
    //// 2. Influenced state - which state is receiving the influence.
    //// 3. The specific step at which the influence is recorded.

    //// Possible values in this vector are:
    //// 1: The state is influencing another state.
    //// -1: The state is being influenced by another state.
    //// 0: Neutral - no influence interaction.
    //// Equivalent to vector<int> blank[numStates][numStates][numSteps+1];
    vector<vector<vector<int>>> blank(numStates, vector<vector<int>>(numStates, vector<int>(numSteps+1, 0)));
    //cout << "Blank vector initialized." << endl;
    //// Choice matrix
    //// choice[current step][current state][influencing step][choice(actual influence or fake influence))]
    //// value is a tuple of (current step, current state, influencing step, choice)

    //tuple<int, int, int> choice[numSteps + 1][numStates][numSteps + 1];
    vector<vector<vector<tuple<int, int, int>>>> choice(numSteps + 1, vector<vector<tuple<int, int, int>>>(numStates, vector<tuple<int, int, int>>(numSteps+1)));

    //cout << "Choice matrix initialized." << endl;


    //// Initial probabilities
    //// init matrix * emission matrix
    //// at time 1 (in a 1-based index system), for all states, there are no influences recorded.
    /// This means that every state at this initial time is in the "0 otherwise" category.
    for (int i = 0; i < numStates; i++){
        m[1][i][1] = E[i][O[1]] * startProb[i];
        I[1][i][1] = blank;
    }

    //// Dynamic programming
    //// Bottom-up approach
    for (int j = 2; j < numSteps + 1; j++){//// current step(from 2 to numSteps). 1 is already calculated(base case).
        for (int r = 0; r < numStates; r++){//// current state
            for (int l = j; l > 0; l--){//// Current influencing step(from j to 1)
                for (int q = 0; q < numStates; q++){//// previous state
                    for (int k = j-1; k > 0; k--){//// previous influencing step(from j-1 to 1)

                        //// Check for affinity rules (if there is an affinity rule, the influence must be next to the influenced state)
                        //// Rule 1: if aff[current state(r)][previous state(q)] == true
                        //// Rule 2: if previous influencing step(k) != current influencing step(l)

                        //// If all rules are satisfied, continue to next iteration
                        //change
                        if (aff[r][q] && k != l+1){
                            continue;
                        }

                        //// Calculate probability
                        //// previous probability and current calculated probability (0 at first)
                        //// Get previous probability (Probability of previous state * transition probability)

                        double prev = m[j-1][q][k] * T[q][r];

                        //// Check if previous probability is 0 (which means it is impossible)
                        if (prev == 0){
                            continue;
                        }
                        double prob = 0;////initialize current probability
                        //// Set current influence matrix to the previous one
                        //// Key part for dynamic programming
                        //// Extract the influence matrix from the previous time step (j-1), given the prior state (q),
                        //// influencing step (k), and current condition (cp). This is the set of all possible influence
                        //// states of the system from the last time step under the provided conditions.
                        vector<vector<vector<int>>> lastInf = I[j-1][q][k];
                        //// Create a working copy of the influence matrix extracted from the previous time step.
                        //// This 'curInf' matrix will serve as the base that will be modified for the current time step (j)
                        //// based on new observations, conditions, and rules. By working with a copy, the original influence
                        //// information from the previous time step is preserved,
                        /// ensuring that the current step's calculations do not affect historical data.

                        vector<vector<vector<int>>> curInf = lastInf;

                        //// Get influencing state
                        //// step, state, influence, choice
                        //// Create a tuple to represent the previous step in terms of step, state, influence, and choice
                        //// next step = (j-1, q, k)

                        tuple<int, int, int> next_step = make_tuple(j-1, q, k);


                        //// Trace path back to influencing step
                        //// Backtrack through the choice matrix to determine the original or earliest influencing step
                        //// that led to the current step (j). This is done by checking the choice matrix at the current step
                        //// (j), current state (r), and current influencing step (l) for the previous step (j-1), previous state (q),
                        //// get<0> means get the first element of the tuple, which is the previous step (j-1)
                        //// l is the current influencing step, so if the previous influencing step is not equal to the current
                        while (get<0>(next_step) > l){
                            //// Get the previous step from the choice matrix
                            //// Choice matrix is a tuple of (current step, current state, influencing step, choice)
                            next_step = choice[get<0>(next_step)][get<1>(next_step)][get<2>(next_step)];
                        }

                        //// Get influencing state
                        //// Extract the state that is causing the influence from the tuple'
                        //// i.e current state (r) = get<1>(next_step)
                        int ingState = get<1>(next_step);

                        //// curinf is the copy of lastinf. Work place.



                        ////add thing
                        //// First if check:
                        //// Rule 1: After influence, the probability is higher than before influence
                        //// Rule 2: Check if influence is valid
                        //// Rule 3: Influence must be at least 4 steps apart
                        //// Rule 4: Check if influenced state is already influenced (-1 is influenced, 1 is influencing, 0 is neutral)
                        if (//// Previous condition is influence, current condition is influence
                                ((prev * N[O[l-1]][O[j-1]]) > prob) && //// After influence, the probability is higher than before influence
                                (validInf(lastInf, ingState, r, l, j)) && ////check if influence is valid
                                (l < j-4) &&//// Influence must be at least 4 steps apart
                                (lastInf[ingState][r][l] != 1)////Check if influenced state is already influenced (-1 is influenced, 1 is influencing, 0 is neutral)
                                ){
                            ////curInf is a 3D vector.
                            ////curInf[current state][influenced state][influencing step]
                            ////curInf[ingState][r][j] = -1 means influenced
                            ////curInf[ingState][r][l] = 1 means influencing
                            curInf[ingState][r][j] = -1;
                            curInf[ingState][r][l] = 1;
                            ////Since the influence is valid, and after influence, the probability is higher than before influence
                            prob = prev * N[O[l-1]][O[j-1]];
                        }

                        //// Fourth if check:
                        //// Rule 1 : cc == 1 (actual)
                        //// Rule 2 : cp == 1 (actual)
                        //// Rule 3 : Current state is not influenceable. Recall Influenceable = {false,true,true,true}
                        //// Rule 4 : current influencing step is equal to previous influencing step, and current state is equal to previous state
                        //// or Rule 4: current influencing step is equal to current step, and current state is not equal to previous state
                        //// Rule 5: Normal HMM Probability (Emmision for observation) is higher than before the influence
                        if ((!influenceable[r]) &&////Current state is not influenceable. Recall Influenceable = {false,true,true,true}
                            ((k == l && r == q) || (l == j && r != q)) &&//// Current influencing step is equal to previous influencing step, and current state is equal to previous state
                            //// or current influencing step is equal to current step, and current state is not equal to previous state
                            ((prev * E[r][O[j]]) > prob)//// Probability (Emmision for observation) is higher than before the influence
                                ){
                            ////normal HMM
                            prob = prev * E[r][O[j]];

                        }
                        //// Recall m is a 4D vector
                        ////m[current step][current state][influencing step][choice(actual influence or fake influence))]
                        ////m[j][r][l][cc] = prob
                        if (prob > m[j][r][l]){
                            ////max probability
                            m[j][r][l] = prob;
                            //cout << "Setting choice[" << j << "][" << r << "][" << l << "] with values (" << j-1 << ", " << q << ", " << k << ")" << endl;  // Debug statement

                            ////choice
                            choice[j][r][l] = make_tuple(j-1, q, k);
                            ////influence matrix
                            I[j][r][l] = curInf;
                        }

                    }
                }
            }
        }


    }

    //// Find highest probability
    ////trace back to find full path
    double maxprob = 0;////Set current highest probability to 0

    //// Creating output strings
    string path[numSteps];//// define a path string with length of numSteps
    string og(numSteps, '.');////initialize string with numSteps * '*', which is unpaired elements
    string influences = og;//// String influences is the same as og

    //// Find highest probabilit end state
    //tuple<int, int, int, int> next_path;
    tuple<int, int, int> next_path;
    int j = numSteps;//// current step j = numSteps
    int longest = 0;////Hidden state path
    //// influences string
    //// For the last step, which is different from other steps
    for (int r = 0; r < numStates; r++){//// current state
        for (int l = j; l >0; l--){//// current influencing step
            //for (int c = 0; c < 2; c++){////current condition (0:fake, 1:actual)
            //if (m[j][r][l][c] > maxprob && endState[r]){////Rule 1: current state is allowed to be end state. Recall: endState is {1,0,0,1}.
            if (m[j][r][l] > maxprob && endState[r]){
                //// L, R are allowed to be end state
                ////Rule 2: current probability is higher than current max probability
                //maxprob = m[j][r][l][c];////update max probability
                //next_path = choice[j][r][l][c];////update next path
                maxprob = m[j][r][l];////update max probability
                next_path = choice[j][r][l];////update next path
                influences = og;////update influences string
                ////r is current state
                ////If R is 1,2,3, which is P, Q, R.
                //// Then influences[l-1] = '(', influences[j-1] = ')'
                if (r == 1 || r == 2 || r == 3){
                    influences[j-1] = ')';
                    influences[l-1] = '(';
                    //// If r is not 0, which is P, Q, R
                    //// Then influences[l-1] = '{', influences[j-1] = '}'
                    ////Wired. Maybe not needed
                } else if (r != 0){
                    influences[j-1] = '}';
                    influences[l-1] = '{';
                }
                ////Hidden state path
                path[j-1] = stateKey[r];
                ////update longest
                longest = stateKey[r].length();
            }
        }

    }
    //iPath = I[][][][]
    // ([current Step] [current state] [influence step] [influence/no influence]) ([influencing state][influenced state][step]) = {1 if influcing, -1 if influenced, 0 otherwise}

    //// traceback to find full path
    ////Recall next_path = choice[j][r][l][c]
    ////next_path = (j-1, q, k, cp)
    for (int j = numSteps-1; j > 0; j--){//// current step
        int r = get<1>(next_path);
        int l = get<2>(next_path);
        //cout << "Accessing choice[" << j << "][" << r << "][" << l << "] for next path." << endl;  // Debug statement

        //int c = get<3>(next_path);
        //if (c == 0){//// if current condition is fake (wired, check later).
        if (r == 1 || r == 2 || r == 3){//// if current state is 1,2,3, which is P, Q, R
            influences[j-1] = ')';
            influences[l-1] = '(';
        } else if (r != 0){//// if current state is not 0, which is P, Q, R. Wired, check later. Maybe not needed
            influences[j-1] = '}';
            influences[l-1] = '{';
            //}

        }
        ////Hidden state path
        path[j-1]= stateKey[r];
        //////update longest
        longest = max(longest, (int)stateKey[r].length());
        //next_path = choice[j][r][l][c];
        next_path = choice[j][r][l];

    }
    vector<string> influenceIndices(influences.length(), "0");  // initialize with "0" strings

// Go through the influence string to find matching brackets and fill their indices
    stack<int> openBrackets;  // To keep track of opening brackets
    for (int i = 0; i < influences.length(); i++) {
        if (influences[i] == '(' || influences[i] == '{') {
            openBrackets.push(i + 1);  // Using 1-based indexing
        } else if (influences[i] == ')' || influences[i] == '}') {
            int openingBracketIndex = openBrackets.top();
            influenceIndices[i] = to_string(openingBracketIndex);
            influenceIndices[openingBracketIndex - 1] = to_string(i + 1); // Using 1-based indexing
            openBrackets.pop();
        }
    }

// Format output

    string actualO = "";
    std::string output;

    // Construct the actualO sequence
    for (auto x : O){
        actualO += observationKey[x];
        for (int j = 0; j < longest - observationKey[x].length() + 1; j++){
            actualO += "";
        }
    }
    output += actualO + "\n";

    // Construct the path sequence
    for (int i = 0; i < numSteps; i++) {
        output += path[i];
        for (int j = 0; j < longest - path[i].length() + 1; j++) {
            output += "";
        }
    }
    output += "\n";

    // Construct the influences sequence
    for (int i = 0; i < influences.length(); i++) {
        output += influences[i];
        for (int j = 0; j < longest; j++) {
            output += "";
        }
    }
    output += "\n";

    return output;


}



int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Please provide a sequence as an argument." << std::endl;
        return 1;
    }
    std::string sequence = argv[1];
    std::string result = predictStructure(sequence);
    std::cout << result << std::endl;
    return 0;
}



