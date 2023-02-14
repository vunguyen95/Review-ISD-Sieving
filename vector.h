#ifndef VECTOR_H
#define VECTOR_H
#include<bits/stdc++.h>
#include<vector>
using namespace std;

/*struct newVec
{
    vector<char> data;
    vector<vector<char>> buckets;

    newVec(): data(0){} 
    newVec(vector<char>& v, int& p);
};*/
//vector<char> mul( vector<char>& b,  vector<vector<char>>& a );
set<vector<int>> samplesGen(set<vector<int>>& res, int& samples, int& length , int& weight);
vector<vector<int>> getComb(vector<int>& v, int n, int& r);
void combination(vector<int>& v, vector<int>& data, int start, int end, int index, int r, vector<vector<int>>& res);

vector<char> represent(vector<int>& v, int& n);
vector<int> bucketRepresent(vector<char>& v);
vector<char> add(vector<char>& u, vector<char>& v);
int weight(vector<char>& v);
vector<vector<char>> parityGen(int& l, int& n);
int rscalar(vector<int>& v1, vector<int>& v2);
vector<int> radd(vector<int>& v1, vector<int>& v2);
vector<int> mul( vector<vector<int>>& a, vector<int>& b);
vector<char> mul1( vector<vector<char>>& a, vector<char>& b);
vector<vector<int>> get(vector<vector<int>>& parity, int& i);
#endif
