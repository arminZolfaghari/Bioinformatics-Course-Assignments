# Bioinformatics Course Coding Assignments

## [HW 1 - Semi-Global Alignment](https://github.com/arminZolfaghari/Bioinformatics-Course-Assignments/tree/main/HW%201%20-%20Semi-Global%20Alignment)
In this assignment, I implemented a ```Semi-Global Alignment``` from scratch. I used ```PAM250 matrix``` for scoring ([PAM250](https://swift.cmbi.umcn.nl/teach/aainfo/pam250.shtml)). 
</br>
The program inputs two protein sequences, calculates the alignment score, and generates two modified sequences after alignment.
#### Sample Input:
```
HEAGAWGHE
PAWHEA
```
#### Sample Output:
```
20
HEAGAWGHE-
---PAW-HEA
```

## [HW 2 - Multiple Sequence Alignment](https://github.com/arminZolfaghari/Bioinformatics-Course-Assignments/tree/main/HW%202%20-%20Multiple%20Sequence%20Alignment)
This assignment is about ```Multiple Sequence Alignment```. [The code](https://github.com/arminZolfaghari/Bioinformatics-Course-Assignments/tree/main/HW%202%20-%20Multiple%20Sequence%20Alignment) consists of two main components. The first part involves the implementation of ```Star Alignment```, a ```Multiple Sequence Alignment``` approach. The second part is a modification step designed to identify and replace blocks in the sequences after alignment.
</br>
The program first prompts the user to input the number of sequences. Subsequently, the sequences are provided as input. The output includes the final ```Multiple Sequence Alignment (MSA)``` score and aligned and modified sequences.
#### Sample Input:
```
4 
TYIMREAQYESAQ
TCIVMREAYE
YIMQEVQQER
WRYIAMREQYES
```
#### Sample Output:
```
51
-TYI-MREAQYESAQ
-TCIVMREA-YE---
--YI-MQEVQQER--
WRYIAMRE-QYES--
```

## [HW 3 - Profile](https://github.com/arminZolfaghari/Bioinformatics-Course-Assignments/tree/main/HW%203%20-%20Profile)
This assignment is about creating a ```profile``` for ```MSA``` and identifying the optimal subsequence of a given sequence.
</br>
The program first prompts the user to input the number of sequences. Subsequently, the sequences are provided as input.
Then, it prompts the user to input a long sequence to find the best subsequence within it. Finally, it prints the best subsequence.
#### Sample Input:
```
4
T-CT
--CT
A-CT
ATCT
ATCCTATATCTTCTCTATACTATCCTTCA
```
#### Sample Output:
```
A-CT
```
