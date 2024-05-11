#include <vector>
#include <set>
#include <iostream>
#include <ctime>
#include <algorithm>


struct WeightPrice {
    WeightPrice(double p=0, const std::vector<int>& w=std::vector<int>()) : price(p), weight(w) {}

    double price = 0;
    std::vector<int> weight;

    WeightPrice operator+(const  WeightPrice& other) const{
        auto tmp = this->weight;
        for (int i = 0; i < tmp.size(); ++i) {
            tmp[i] += other.weight[i];
        }
        return {this->price + other.price, tmp};
    }

    WeightPrice operator-(const  WeightPrice& other) const{
        auto tmp = this->weight;
        for (int i = 0; i < tmp.size(); ++i) {
            tmp[i] -= other.weight[i];
        }
        return {this->price - other.price, tmp};
    }

    bool operator<(const  WeightPrice& other) const{
            return (this->price < other.price) || (this->price == other.price) && (this->weight[0] < other.weight[0]);
        }

};


class NDimensionalBranch {
public:
    void Init(const std::vector<int>& max_w, const std::vector<int>& prices, const std::vector<std::vector<int>>& weights) {
        max_weight_ = max_w;
        input_.resize(prices.size());
        for (int i = 0; i < input_.size(); ++i) {
            input_[i] = {prices[i], weights[i]};
        }
        cur_price_weight_ = {0, std::vector<int>(max_w.size(), 0)};
    }
    
    bool CheckWeightSum(const std::vector<int>& first, const std::vector<int>& second) {
        for (int i = 0; i < second.size(); ++i) {
            if (first[i] + second[i] > max_weight_[i]) {
                return false;
            }
        }
        return true;
    }

    int Solve() {
        Dfs(input_.size() - 1, true);
        Dfs(input_.size() - 1, false);
        return best_price_;
    }

    void Dfs(int node, bool pick) {
        // static int count = 0;
        // ++count;
        // if (count % 1000 == 0) {
        //     std::cout << count << std::endl;
        // }
        ++n;
        if (node < 0) {
            return;
        }
        if (pick && !CheckWeightSum(cur_price_weight_.weight, input_[node].weight)) {
            return;
        }
        if (pick) {
            cur_price_weight_ = cur_price_weight_ + input_[node];
            cur_picked_.push_back(node);
        }
        if (cur_price_weight_.price > best_price_) {
            best_price_ = cur_price_weight_.price;
            best_pick_ = cur_picked_;
        }   
        // auto tmp = SolveLagrange(node - 1);
        // if (tmp + cur_price_weight_.price > best_price_) {
        //     // std::cout << "NO BLOCK" << std::endl;
        Dfs(node - 1, true);
        Dfs(node - 1, false);
        // }
        // } else {
        //     std::cout << "blocked" << std::endl;
        // }
        if (pick) {
            cur_price_weight_ = cur_price_weight_ - input_[node];
            cur_picked_.pop_back();
        }
    }
int n = 0;
private:
    std::vector<WeightPrice> input_;
    std::vector<int> cur_picked_;
    std::vector<int> max_weight_;
    WeightPrice cur_price_weight_;
    std::vector<int> best_pick_;
    int best_price_ = -1;
};
