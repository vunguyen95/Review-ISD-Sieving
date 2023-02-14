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
    /*for (int j = 0; j< weight; j++)
    {
        auto position = distrib(gen);
        while(find(temp.begin(),temp.end(),position)!= temp.end())
    	{
            position = distrib(gen);
        }
    temp.emplace_back(position);
    }*/
	sort(temp.begin(), temp.end());
	return temp;
}
//TEST CLEAR 
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
void sort_samples(Set& s, Map& m) // put each vector according to its first coordinate
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
void combine(Map& m, Set& s, int p, int two_p, Vec pparity, vector<int> psyndrome, int round, int nbrsamples, int length)
{
    int p_minus = p-1;
    vector<int> zero(round);
    int key = 1 ;
    //float count{0};
    //float duplicate0{0};
    //float duplicate1{0};
    //float count_total{0};
    //Set s0 = s; // Check
    //Set s00;
    //check_Map1 m00;
    //Set s1; // Combine
    //Set s11;
    //check_Map1 m11;
    //cout << count << ": " << ": " << duplicate0 << ": " << duplicate1 << endl;
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
                    //cout << (*begin).first << endl;
                    //cout << (*begin).second << endl;
                    //std::cout << (*it) << endl;
                    std::cerr << e.what() << "here1" << '\n';
                }
                
                // A sample put it to a label many times here, fixed it but still dont know why.
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
			            //if(s.size() >= nbrsamples){break;} //|| count >= nbrsamples)
                        for(int j = i+1; j < (*it).second.size(); j++)
                        {
				            //if(s.size() >=nbrsamples){break;} //|| count >= nbrsamples
                            auto sum = radd((*it).second[i], (*it).second[j]);
                            if(sum.size() == two_p)
                            {

                                auto check = mul(pparity,sum);
                                if(check == psyndrome || check == zero)
                                {	//count_total++;
                                    if(s.find(sum) == s.end()){
                                        s.insert(sum);
                                    }
                                    /*if(s0.find(sum) == s0.end() && s1.find(sum) == s1.end()){
                                        count++;
                                        s1.insert(sum);
                                    }
                                    else if(s0.find(sum) != s0.end() && s1.find(sum) == s1.end()){
                                        duplicate0++;
                                        s00.insert(sum);
                                        auto pround = round-1;
                                        vector<int> ppsyndrome;
                                        ppsyndrome.assign(psyndrome.begin(), psyndrome.begin() + pround);
                                        vector<int> ppzero(round-1);
                                        auto check0 = (*it).second[i];
                                        auto ppparity = get(pparity, pround);
                                        if(mul(ppparity, check0) == ppsyndrome ||mul(ppparity, check0) == ppzero){
                                            check0.insert(check0.end(), 1);
                                        }
                                        else{
                                            check0.insert(check0.end(), 0);     
                                        }
                                        auto check1 = (*it).second[j];
                                        if(mul(ppparity, check1) == ppsyndrome ||mul(ppparity, check1) == ppzero){
                                            check1.insert(check1.end(), 1);
                                        }
                                        else{
                                            check1.insert(check1.end(), 0);     
                                        }
                                        if(m00.find(sum) == m00.end()){
                                            vector<vector<int>> pairs;
                                            pairs.emplace_back(check0);
                                            pairs.emplace_back(check1);
                                            m00.insert(make_pair(sum, pairs));
                                        }
                                         else{
                                            auto it2 = m00.find(sum);
                                            (*it2).second.emplace_back(check0);
                                            (*it2).second.emplace_back(check1);
                                        }
                                    }
                                    else if(s1.find(sum) != s1.end() && s0.find(sum) == s0.end()){
                                        duplicate1++;
                                        s11.insert(sum);

                                        if(m11.find(sum) == m11.end()){

                                            vector<vector<int>> Pairs;
                                            auto check = (*it).second[i];
                                            check.insert(check.end(),(*it).first.begin(), (*it).first.end());
                                            Pairs.emplace_back(check);
                                            Pairs.emplace_back((*it).second[j]);
                                            m11.insert(make_pair(sum, Pairs));
                                        }
                                        else{
                                            auto it1 = m11.find(sum);
                                            auto check = (*it).second[i];
                                            check.insert(check.end(),(*it).first.begin(), (*it).first.end());
                                            (*it1).second.emplace_back(check);
                                            (*it1).second.emplace_back((*it).second[j]);
                                        }
                                    }*/
                                    
                                    
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
    /*cout << "total collision: " << count_total << " log: " << log2(count_total) << endl;
    cout << "Unique Collision Count: " << log2(count) << endl;
    cout << "Duplicate from preliminary check: " << log2(duplicate0) <<" : " <<duplicate0 << endl;
    //cout << "Number of duplicated vectors from check: " << s00.size() << endl;
    cout << "Duplicate from combination: " << log2(duplicate1) <<" : " <<duplicate1 <<endl;
    //cout << "Number of duplicated vectors from combination: " << s11.size() << endl;
    cout << "Total duplicates:" << duplicate0 + duplicate1 << endl;
    cout << " Check Total - unique - duplicate: " << count_total - count - duplicate0 - duplicate1 << endl;
    if(duplicate0 + duplicate1 > 0){
    cout << "Ratio Duplicate/M: " << static_cast<float>((duplicate0 + duplicate1)/nbrsamples) << endl;
    }*/
    /*vector<vector<int>> max0;
    vector<int> vec0;
    vector<int> check_syndrome;
    vector<vector<int>> max1;
    vector<int> vec1;
    for(auto it = m00.begin(); it!= m00.end();it++){
        if((*it).second.size() > max0.size()){
            max0 = (*it).second;
            vec0 = (*it).first;
            check_syndrome = mul(pparity, vec0);
        }
    }
    for(auto it = m11.begin(); it!= m11.end(); it++){
        if((*it).second.size() > max1.size()){
            max1 = (*it).second;
            vec1 = (*it).first;
        }
    }
    cout << "Most duplicated vector from check: " << vec0;
    cout << max0 << endl;
    cout << "Check syndrome:" << check_syndrome; 
    cout << "Most duplicated vector from combination: " << vec1;
    cout << max1 << endl;*/
}

int main()
{
   //Parameters:
   int length = 1000+ 50 ;// k+ l
   int nbrParities = 50;
   int p = 2;
   int two_p = 4;
   int nbrSamples = pow(2,15.406);
   int count_success = 0;
   int attemps = 100;
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
	check(myMap, partialParity, partialSyndrome, mySet,round);
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

