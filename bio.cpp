#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <iterator>
#include "bio.h"

using std::vector, std::istringstream, std::string, std::isupper, std::swap, std::find, std::begin, std::end, std::distance;

//function that returns the RNA sequence but not the reverse complement. 
//used in the function GetReadingFramesAsCodons. see line 87.
string RNAOriginal(string input){
  
  for(char &i: input){
    if(i == 'T'){
      i = 'U';
    }
  }

  return input;
}

bool IsValidDNASequence(const std::string & input){
  for(const char &i: input){
    //referenced the hw 8 solution for std::isupper
    if(!isupper(i)){
      return false;
    }
    else if(i != 'A' && i != 'T' && i != 'C' && i != 'G'){
      return false;
    }
  }
  return true;
}

void GetReverseComplementSequence(const std::string & input,  std::string * const output){
  istringstream iss(input);
  char ch = ' ';
  string RNA_sequence = "";

  while(iss >> ch){
    if (ch == 'A' && *output != "rna"){
      RNA_sequence.insert(0, 1, 'T');
    }
    else if (ch == 'T'){
      RNA_sequence.insert(0, 1, 'A');
    }
    else if (ch == 'C'){
      RNA_sequence.insert(0, 1, 'G');
    }
    else if (ch == 'G'){
      RNA_sequence.insert(0, 1, 'C');
    }
    else{
      RNA_sequence.insert(0, 1, 'U');
    }
  }
  
  *output = RNA_sequence;


}

std::string GetRNATranscript(const std::string & input){
  //assigned output to "rna" so it could be recognized as a 
  //RNA sequence in GetReverseComplementSequence. See line 43.
  string output = "rna";

  GetReverseComplementSequence(input, &output);

  return output;
}

std::vector<std::vector<std::string>> GetReadingFramesAsCodons(const std::string & input){
  vector<vector<string>> codons_vectors(6);
  string output = "";
  string sequence = "";
  string temp = "";
  int count = 0;

  sequence = GetRNATranscript(input);
  
  //the for loop will iterate through each vector<string> inside codons
  for (auto &codons : codons_vectors){
    if(count == 3){ //once loop iterates 3 times the sequence is overwritten with the no reversed RNA sequence
      sequence = RNAOriginal(input);
    }
    //the for loop will insert a string with a size = 3 into a vector's indexes.
    for(char i: sequence){
      temp.push_back(i);
      if(static_cast<int>(temp.size()) == 3){
        codons.push_back(temp);
        temp = "";
      }
    }
    temp = "";
    sequence.erase(0, 1);
    ++count;
  }
  
  
  return codons_vectors;
}

std::string Translate(const std::vector<std::string> & codon_sequence){
  string temp = "";
  vector<string> codon_list = {"GCU", "GCC", "GCA", "GCG", "CGU", "CGC", "CGA", "CGG", "AGA", "AGG",
  "AAU", "AAC", "GAU", "GAC", "UGU", "UGC", "CAA", "CAG", "GAA", "GAG",
  "GGU", "GGC", "GGA", "GGG", "CAU", "CAC", "AUU", "AUC", "AUA", "UUA",
  "UUG", "CUU", "CUC", "CUA", "CUG", "AAA", "AAG", "AUG", "UUU", "UUC",
  "CCU", "CCC", "CCA", "CCG", "UCU", "UCC", "UCA", "UCG", "AGU", "AGC",
  "ACU", "ACC", "ACA", "ACG", "UGG", "UAU", "UAC", "GUU", "GUC", "GUA", 
  "GUG", "UAG", "UGA", "UAA"};
  vector<string> amino_acids = {"A", "A", "A", "A", "R", "R", "R", "R", "R", "R", "N", "N", "D", "D",
  "C", "C", "Q", "Q", "E", "E", "G", "G", "G", "G", "H", "H", "I", "I",
  "I", "L", "L", "L", "L", "L", "L", "K", "K", "M", "F", "F", "P", "P",
  "P", "P", "S", "S", "S", "S", "S", "S", "T", "T", "T", "T", "W", "Y",
  "Y", "V", "V", "V", "V", "*", "*", "*"};
  
  //I looked over these two sources to help with the for loop:
  //https://piazza.com/class/k4wvjqt9b2x72k?cid=786
  //https://thispointer.com/c-how-to-find-an-element-in-vector-and-get-its-index/
  for(const string &codon : codon_sequence){
    auto pos = find(codon_list.begin(), codon_list.end(), codon);
    auto dist = distance(begin(codon_list), pos);
    temp.append(amino_acids[dist]);

  }
  
  return temp;
}

std::string GetLongestOpenReadingFrame(const std::string & DNA_sequence){
  string translated_sequence = "";
  string longest_frame = "";
  string temp = "";
  vector<vector<string>> codons_vectors = GetReadingFramesAsCodons(DNA_sequence);
  
  for(vector<string> codons: codons_vectors){
    translated_sequence.append(Translate(codons)); //translate the codons to amino acids and append it to translated_sequence
    int start = static_cast<int>(translated_sequence.find('M'));
    int end = static_cast<int>(translated_sequence.find('*', start)); //search for '*' from where we found 'M' to the end
    
    //while loops ends once there is no 'M' or '*' left in the sequence
    while(start != -1 && end != -1){
      temp = translated_sequence.substr(start, end - start + 1); //potential longest open reading frame
      if(static_cast<int>(temp.size()) > static_cast<int>(longest_frame.size())){
        longest_frame = temp;
      }
      translated_sequence.erase(0, end + 1); //erase the substring from 'M' to '*'
      start = static_cast<int>(translated_sequence.find('M'));
      end = static_cast<int>(translated_sequence.find('*', start));
    }
    
    translated_sequence = "";
  }

  return longest_frame;
}