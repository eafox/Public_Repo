# Miscellaneous Work

#### Background
This is a challenge I completed in 2018 that my interviewer has graciosuly allowed me to share. 

## Prompt:
Created on Wed Jul 18 09:07:53 2018

@author: ofinkler

Xyz's biggest client needs us to help him increase his production
efficiency by predicting the end state of a given ore sample.
The production procedure starts as raw ore, then during processing, the ore 
begins randomly changing between forms, eventually reaching a stable form 
(also reffered to as a terminal state). 

The probability for each ore structure to transform is fixed.
That means that each time the ore is in one state it has the same probabilities 
of changing to the next state. Note, that the ore might also has a certain 
probability of staying at the same state. 
The observed state transitions are recorded in a matrix. 

Your assignment is to write a script contianing a function solution(input) that 
takes an array of array of nonnegatvie integers, representing the number of 
times each state has gone to the next state.
The function should return an array of integers for each stable state (terminal 
state), giving the exact probabilities of each terminal state, represented as 
the numerator for each state and then the denominator for all of them at the end 
and in simplest form (see example below). 

The matrix is at most 10 by 10. 

The problem assumes that no matter which state the ore is in, there
is a path from that state to a terminal state. 
This means that all processing will eventually end in a stable state.

The ore always starts in state 0. 
The denominator will fit within a signed 32-bit integer during the calculation, 
as long as the fraction is simplified regularly.

For example, consider the matrix m:\
[ \
 [0,1,0,0,0,1], # s0, the inital state, goes to s1 and s5 with equal probability\
 [4,0,0,3,2,0], # s1 can become s0, s3, or s4, but with different probabilities\
 [0,0,0,0,0,0], # s2 is terminal, and unreachable (never observed in practice)\
 [0,0,0,0,0,0], # s3 is terminal\
 [0,0,0,0,0,0], # s4 is terminal\
 [0,0,0,0,0,0]  # s5 is terminal\
]

So, we consider different paths to terminal states, such as:
s0 -> s1 -> s3\
s0 -> s1 -> s0 -> s1 -> s0 -> s1 -> s4\
s0 -> s1 -> s0 -> s5\

Tracing the probabilities of each, we find that \
s2 has probability 0\
s3 has probability 3/14\
s4 has probability 1/7\
s5 has probability 9/14 

So, putting that together, and making a common denominator, gives an answer in 
the form of \
[s2.numerator, s3.numerator, s4.numerator, s5.numerator, denominator] \
which is \
[0, 3, 2, 9, 14]

Test cases\
Inputs: (int) m = [[0, 2, 1, 0, 0], [0, 0, 0, 3, 4], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]]

Output: (int list) [7, 6, 8, 21]

Inputs: (int) m = [[0, 1, 0, 0, 0, 1], [4, 0, 0, 3, 2, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0]]

Output: (int list) [0, 3, 2, 9, 14]

## My Response:
_see interview_challenge_2018-EmmaFox.py for code_

The script contains a function (and its dependent) functions to solve for the final distribution between terminal states. The script can be run \
A) On its own from the terminal (ex: python3 UC_challenge_2018-EmmaFox.py) or python command-line which will print out the answers to the test cases. \
B) From the terminal with an input matrix. Ex: \
_python3 UC_challenge_2018-EmmaFox.py "[[0, 2, 1, 0, 0], [0, 0, 0, 3, 4], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]]"_\
This will print out the final distribution for that particular input matrix. \
C) From within python either using the step-by-step code at the bottom of the file or using the full pipeline function "main()" with a custom input matrix. 

I chose to solve the problem by taking advantage of the fact that the system behaves as an absorbent Markov chain. I had some experience with transition matrices because they're used in ecology to calculate the stable age-structure of a population. Though my previous experience was mostly with eigen-analysis, I was able to expand upon this with further linear algebra research in order to determine which type of problem this situation fit into. 

For cases with no loops/cycles, like test case 1, the final distribution could have been calculated by simply multiplying the transition probabilities along each path and summing across paths to that particular terminal state. However, this would require a test to check whether loops existed and extra modification to handle cycles like the one between state 0 and state 1 in test case 2. On the other hand, working with matrices for absorbent Markov chains can handle both cases without extra modification and still return the answers in a comparable amount of time. 

At the moment, the script returns a list of the distribution between all terminal states as floats. It is not possible to complete the matrix operations needed for the formulas (particularly the inverse of a large matrix) using fraction objects rather than floats in python at this time, according to my research. The rounding that occurs with operations on floats meant that converting from the final distribution floats to fractions directly gave denominators larger than 10^9 in most cases.  For cases without loops, like test case 1, the denominator could be calculated as the product of the denominator from the two non-terminal stages. For cases with loops, the denominator could be calculated as the product of the denominator of the two non-terminal stages minus the numerator of the fraction (with the same denominator) that was likely to remain in the cycle after 1 run. For example, in test case 2, the product of the two non-terminal denominators was 18 (2*9). After 1 iteration, the likely distribution is 9/18 in terminal state 5, 3/18 in terminal state 3, and 2/18 in terminal state 4 with 4/18 having transitioned back from state 1 to state 0. So, the denominator would be 14 (18-4). While I noted this relationship, I ran out of time to test whether this would hold true with more complex systems (i.e. more than 2 non-terminal states, presence of multiple loops) and devise a function to convert the output. However, I did check my float answers against the test answers provided and they did match.