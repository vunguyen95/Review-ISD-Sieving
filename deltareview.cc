#include"misc.h"
#include"vector.h"
#include<algorithm>
#include<random>
#include<cstdlib>
#include<unordered_set>
#include<map>
#include<unordered_map>
#include<chrono>
using namespace std;
//This is to show that without the checking part, the number of combination is consistent with (14), without $\delta$
// Modified: comment out checking part
// Remember that we have to increase M by a factor of 2 (because we did not keep any vector from the previous iteration.)
template < typename SEQUENCE > struct seq_hash //Hash function for vectors
{
	long long int operator() (const SEQUENCE& seq) const
	{
		std::size_t hash = 0;
		//boost::hash_range(hash, seq.begin(), seq.end());
		for (auto &i : seq) { hash ^= i + 0x9e3779b9 + (hash << 6) + (hash >> 2); }
			return hash;
	}
};
using namespace std::chrono;
typedef unordered_map<int, vector<vector<int>>> Map;
using temp_Map = std::unordered_map< vector<int>, vector<vector<int>>, seq_hash<vector<int>> >; 
typedef vector<vector<int>> Vec;
using check_Map = std::unordered_map< vector<int>, int, seq_hash<vector<int>> >;
using check_Map1 = std::unordered_map< vector<int>, vector<vector<int>>, seq_hash<vector<int>> >;
using Set = std::unordered_set< vector<int>, seq_hash<vector<int>> >;

/* Generating random weight-p vectors */
/*
length: length of vectors
weight: weight of vectors
*/
vector<int> vector_gen(int& length, int& weight)
{
	random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<int> distrib(1,length);
    vector<int> temp(weight);
    set<int> s;
    int i = 0;
    while(i < weight){
        int position = distrib(gen);
        if (s.count(position) == 0){
            temp[i++] = position;
            s.insert(position);
        }
    }
	sort(temp.begin(), temp.end());
	return temp;
}
//TEST CLEAR 

/* Generating list of vectors */ 
/*
res: destination container. Here we use set<vector<int>>

*/
Set sample_gen(Set& res, int& samples, int& length , int& weight)
{
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<int> distrib(1,length);
    while(res.size()<= samples)
    {
        auto temp = vector_gen(length, weight);
        res.insert(temp);
    }
    return res;
}

// TEST CLEAR 
/* Move a vector to a `bucket'. Here the key is an integer corresponding to p' = 1. */
void move(Map& m, vector<int> sample, int& key)
{
    auto find = m.find(key);
    if(find == m.end())
    {
        Vec bucket; 
        bucket.emplace_back(sample);
        m.insert(pair<int, Vec>(key, bucket));
       
    }
    else
    {
        if(std::find((*find).second.begin(), (*find).second.end(), sample) == (*find).second.end()){
                (*find).second.emplace_back(sample);
            }
    }
}
// TEST CLEAR 

/* After generating the sample, we put each vector to `buckets' according to its first coordinate*/
void sort_samples(Set& s, Map& m) // 
{
    for(auto sample : s)
    {
        int key = sample[0];
        auto find = m.find(key);
        if(find == m.end())
        {
            Vec bucket;
            bucket.emplace_back(sample);
            m.insert(pair<int, Vec>(key, bucket));
        }
        else
        {
            if(std::find((*find).second.begin(), (*find).second.end(), sample) == (*find).second.end()){
                (*find).second.emplace_back(sample);
            }
            
        }
    }
    s.clear();

}
//TEST CLEAR 

/* Finding the next bucket of a sample, then move it*/
void move_samples(Map& m, int& key, int p)
{
    auto it = m.find(key);
    for(auto sample : (*it).second)
    {
        auto next = find(sample.begin(), sample.end(), key)+1;
        if(*next <= sample[p]) //stop
        {   
            auto next_key = *next;
            move(m, sample, next_key);
        }
    }
}

/* Checking if a vector satisfies the partial syndrome condition*/
void check(Map& m, Vec& partialParity, vector<int>& partialSyndrome, Set& s, int& round) //check
{
	//vector<int> zero(round);
	vector<int> zero(partialSyndrome.size());
	for(auto it = m.begin(); it!= m.end(); it++)
	{
		for(int i = 0; i < (*it).second.size(); i++)
		{
			auto syndrome = mul(partialParity, (*it).second[i]);
			if(syndrome == partialSyndrome || syndrome == zero)
			{
				s.insert((*it).second[i]);
			}
		}
	}
}
//

/* At bucket label p', we use a hash table indexed by p'' = p - p' to detect collisions*/
void combine(Map& m, Set& s, int p, int two_p, Vec pparity, vector<int> psyndrome, int round, int nbrsamples, int length)
{
    int p_minus = p-1; 
    vector<int> zero(round);
    int key = 1 ;
    
    while(key<= length && m.size()> 0 && s.size() < nbrsamples)
    {    
        auto begin = m.find(key);
        if(begin!= m.end())
        {
            temp_Map table;
            for(auto sample: (*begin).second)
            {
                auto it = find(sample.begin(),sample.end(),key)+1;
                if(it!= sample.end())
                {
                    vector<int> remain;
                try
                {
                    remain.assign(it, sample.end());
                    
                }
                catch(const std::exception& e)
                {
                  
                    std::cerr << e.what() << "here1" << '\n';
                }
                
                auto list = getComb(remain, remain.size(), p_minus);
                for (auto label: list)
                {
                    auto find = table.find(label);
                    if(find == table.end())
                    {
                        Vec new_bucket;
                        new_bucket.emplace_back(sample);
                        table.insert(make_pair(label, new_bucket));
                    }
                    else
                    {
                        if(std::find((*find).second.begin(), (*find).second.end(), sample) == (*find).second.end()){
                            (*find).second.emplace_back(sample);
                        }
                    }
                }
                }
                else
                {
                    cout << (*begin).first << endl;
                    cout << sample;
                }

            }
            for(auto it = table.begin(); it!= table.end();it++)
            {
                if((*it).second.size() > 1)
                {
                    for (int i = 0; i < (*it).second.size()-1; i++)
		    {
			            if(s.size() >= nbrsamples){break;} //|| count >= nbrsamples)
                        for(int j = i+1; j < (*it).second.size(); j++)
                        {
				            if(s.size() >=nbrsamples){break;} //|| count >= nbrsamples
                            auto sum = radd((*it).second[i], (*it).second[j]);
                            if(sum.size() == two_p)
                            {

                                auto check = mul(pparity,sum);
                                if(check == psyndrome || check == zero)
                                {	//count_total++;
                                    if(s.find(sum) == s.end()){
                                        s.insert(sum);
                                    }
                                    
                                }
                            }
                        }
                    }
                }
            }
            table.clear();
            move_samples(m, key, p);
            m.erase(begin);
            key++;
    
        }
        else
        {
            key++;
        }

    }
    
}

int main()
{
   //Parameters:
   int length = 300 + 28 ;// k+ l
   int nbrParities = 28;
   int p = 3;
   int two_p = 6;
   int nbrSamples = pow(2,15.5); //No delta here, M = 4/q because of no checking.
   int count_success = 0;
   int attemps = 1;
   for(int run = 0; run < attemps; run++)
   {
	Map myMap;
	Set mySet;
    cout << length<< endl;
    //Parity, random test error e and syndrome:
   auto h0 = parityGen(nbrParities,length);

   vector<vector<int>> H0(h0.size());
		for (int i = 0; i < H0.size(); i++)
		{
			H0[i] = bucketRepresent(h0[i]);
		}


   auto testError = vector_gen(length, two_p);

   cout << testError;
   auto syndrome = mul(H0, testError); 
   auto check1 = represent(testError,length);
   //cout << "Parity check H0: " <<endl << h0;
   //cout << "Bucket form: " << endl << H0;
   //cout << "Test error: " << endl << testError;
   cout << "Syndrome:" << endl << syndrome;
   cout << mul1(h0,check1);

   //ALGORITHM
   sample_gen(mySet, nbrSamples, length, two_p);
   sort_samples(mySet,myMap);

   for(int round = 1; round <= nbrParities; round++)
   {
	cout << "ROUND" <<round <<endl;
   	sort_samples(mySet,myMap);
   	auto partialParity = get(H0, round);
   	vector<int> partialSyndrome;
   	partialSyndrome.assign(syndrome.begin(), syndrome.begin() + round);
    //cout << partialSyndrome;
	//check(myMap, partialParity, partialSyndrome, mySet,round);
	//cout << "Preliminary check: "<<log2(mySet.size()) <<endl;

   	combine(myMap, mySet,p, two_p, partialParity, partialSyndrome, round, nbrSamples, length);
	cout << "List size after " << round << "iteration: " <<log2(mySet.size()) << endl;
   }
   auto it = mySet.find(testError);
   if(it!= mySet.end())
   {
	   cout << "FOUND SOLUTION" << endl;
	   //cout << testError;
	   //cout << *it;
       	count_success++;
   }
   else
   {cout << "SOLUTION NOT FOUND" << endl;}
   }
   cout <<"Count success: "<< count_success << endl;
}

