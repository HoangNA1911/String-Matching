#include <iostream>
#include <string>
#include <fstream>
#include <cstring>
#include <vector>
using namespace std;

#define d 256
struct wordIF
{
    string name;
    int x;
    int y;
    string dir;
};
void print(char **matrix, int w, int h);
void readFile(char **&matrix, int &width, int &height, vector<string> &word);
void readFile(char **&matrix, int &width, int &height, vector<string> &word)
{
    ifstream fin;
    fin.open("input.txt", ios::in);
    string s;
    getline(fin, s, ' ');
    width = stoi(s);
    getline(fin, s);
    height = stoi(s);
    matrix = new char *[height];
    for (int i = 0; i < height; i++)
    {
        matrix[i] = new char[width];
    }
    for (int i = 0; i < height; i++)
    {
        int k = 0;
        // char ch;
        for (int j = 0; j < width - 1; j++)
        {
            getline(fin, s, ' ');
            // strcpy(&ch, s.c_str());
            matrix[i][k] = s[0];
            k++;
        }
        getline(fin, s);
        // strcpy(&ch, s.c_str());
        matrix[i][k] = s[0];
    }
    getline(fin, s);
    while (!fin.eof())
    {
        getline(fin, s);
        if (s == "#")
            break;
        word.push_back(s);
    }
    fin.close();
}
void print(char **matrix, int w, int h)
{
    for (int i = 0; i < h; i++)
    {
        for (int j = 0; j < w; j++)
        {
            cout << matrix[i][j] << " ";
        }
        cout << endl;
    }
}
vector<wordIF> brute_force(char **matrix, int width, int height, vector<string> &word)
{
    vector<wordIF> result;
    wordIF dummy;
    int k = 0;
    while (k < word.size())
    {
        string s = word[k];
        int M = s.length();

        for (int i = 0; i < height; i++)
        {
            for (int j = 0; j <= width - M; j++)
            {
                int h;
                for (h = 0; h < M; h++)
                {
                    if (matrix[i][j + h] != s[h])
                        break;
                }

                if (h == M)
                {
                    dummy.x = j;
                    dummy.y = i;
                    dummy.dir = "LR";
                    dummy.name = word[k];
                    result.push_back(dummy);
                }
            }
        }
        for (int j = 0; j < width; j++)
        {
            for (int i = 0; i <= height - M; i++)
            {
                int h;
                for (h = 0; h < M; h++)
                {
                    if (matrix[i + h][j] != s[h])
                        break;
                }
                if (h == M)
                {
                    dummy.x = j;
                    dummy.y = i;
                    dummy.dir = "TD";
                    dummy.name = word[k];
                    result.push_back(dummy);
                }
            }
        }
        k++;
    }
    return result;
}



// ham Rabin karp

vector<wordIF> rabin_Karp(char **matrix, int width, int height, vector<string>& word, int q)
{
    vector<wordIF> result;
    wordIF dummy;
    int k = 0;
    while (k < word.size())
    {
        string s = word[k];
        int M = s.length();
        for (int i = 0; i < height; i++)
        {
           unsigned long long h = 1;
            for (int p = 0; p < M - 1; p++)
                h = (h * d) % q;
            // first hash
            long long a = 0, b = 0;
            for (int p = 0; p < M; p++)
            {
                a = (a * d + s[p]) % q;
                b = (b * d + matrix[i][p]) % q;
            }
            for (int j = 0; j <= width - M; j++)
            {
                int p;
                if (a == b)
                {
                    for (p = 0; p < M; p++)
                    {
                        if (matrix[i][j + p] != s[p])
                            break;
                    }
                    if (p == M)
                    {
                        dummy.x = j;
                        dummy.y = i;
                        dummy.dir = "LR";
                        dummy.name = word[k];
                        result.push_back(dummy);
                    }
                }
                if (j < width - M)
                {
                    b = (d * (b - matrix[i][j] * h) + matrix[i][j + M]) % q;
                    if(b<0) b=b+q;
                }
            }
        }
        for (int j = 0; j < width; j++)
        {
           unsigned long long h = 1;
            for (int p = 0; p < M - 1; p++)
                h = (h * d) % q;
            // first hash
            long long a = 0, b = 0;
            for (int p = 0; p < M; p++)
            {
                a = (a * d + s[p]) % q;
                b = (b * d + matrix[p][j]) % q;
            }
            for (int i = 0; i <= height - M; i++)
            {
                int p;
                if (a == b)
                {
                    for (p = 0; p < M; p++)
                    {
                        if (matrix[i+p][j] != s[p])
                            break;
                    }
                    if (p == M)
                    {
                        dummy.x = j;
                        dummy.y = i;
                        dummy.dir = "TD";
                        dummy.name = word[k];
                        result.push_back(dummy);
                    }
                }
                if (j < height - M)
                {
                    b = (d * (b - matrix[i][j] * h) + matrix[i][j + M]) % q;
                }
            }
        }
        k++;
    }
    return result;
}

// ham KMP 
void LBS(string word,int lps[]){
    lps[0]=0;
    int idx=0;
    int i=1;
    while(i<word.size()){
        if(word[i]==word[idx]){
            idx++;
            lps[i]=idx;
            i++;
        }
        else{
            if(idx>0){
                idx=lps[idx-1];
            }
            else{
                lps[i]=0;
                i++;
            }
        }
    }
}
vector<wordIF> KMD(char **matrix, int width, int height, vector<string> &word){
    vector<wordIF> result;
    wordIF dummy;
    int k=0;
    while(k<word.size()){
        string s=word[k];
        int M=s.size();

        for(int i=0;i<height;i++){
            int lps[M];
            LBS(s,lps);
            int h=0;
            int j=0;
            while(j<width){
                if(matrix[i][h]==s[j]){
                    h++;
                    j++;
                }
                if(j==s.size()){
                    dummy.x = j;
                    dummy.y = i;
                    dummy.dir = "TD";
                    dummy.name = word[k];
                    result.push_back(dummy);
                    j=lps[j-1];
                }
                else if(h<width&& matrix[i][h]!=s[j]){
                    if(j!=0) j=lps[j-1];
                    else{
                        i=i+1;
                    }
                }
            }
        }
        

        k++;
    }
    return result;
}

void printWPD(vector<wordIF> result)
{
    cout << result.size() << endl;
    for (int i = 0; i < result.size(); i++)
    {
        cout << result[i].name << " (" << result[i].x << "," << result[i].y << ") " << result[i].dir << endl;
    }
}
int main()
{
    // doc file nguoi dung luu vao mang 2 chieu cap phat dong
    // bao gom: doc so hang va so cot (cap phat dong cho mang 2 chieu va de biet can doc bao nhieu hang trong file de luu vao mang 2 chieu).
    //          doc cac ky tu sau do luu vao mang (luu vao 1 mang 2 chieu kieu char)
    //          doc cac tu can tim kiem trong mang. (luu vao 1 mang 1 chieu string)
    char **matrix;
    int width, height;
    vector<string> word;
    readFile(matrix, width, height, word);
    print(matrix, width, height);
    // tim kiem tren mang, tim kiem theo 3 thuat toan Brute-force, Rabin-Karp, Knuth-Morris-Pratt
    vector<wordIF> checkMatching = brute_force(matrix, width, height, word);
    printWPD(checkMatching);
    int q = INT_MAX;
    checkMatching= KMD(matrix, width, height, word);
    printWPD(checkMatching);

    for (int i = 0; i < height; i++)
    {
        delete[] matrix[i];
    }
    delete[] matrix;
    return 0;
}
