#include <vector>
#include <limits>
#include <set>
#include <iostream>


class NDimensionalKnapsack {
public:
    struct WeightPrice {
        WeightPrice(int p=0, const std::vector<int>& w=std::vector<int>()) : price(p), weight(w) {}

        int price = 0;
        std::vector<int> weight;

        static bool VectorLess(const std::vector<int>& one, const std::vector<int>& other) {
            if (one == other) {
                return false;
            }
            for (int i = 0; i < one.size(); ++i) {
                if (one[i] > other[i]) {
                    return false;
                }
            }
            return true;
        }

        bool operator<(const  WeightPrice& other) const{
            return (this->price < other.price) || (this->price == other.price) && VectorLess(this->weight, other.weight);
        }

        WeightPrice operator+(const  WeightPrice& other) const{
            auto tmp = this->weight;
            for (int i = 0; i < tmp.size(); ++i) {
                tmp[i] += other.weight[i];
            }
            return {this->price + other.price, tmp};
        }

    };

    void Init(const std::vector<int>& max_w, const std::vector<int>& prices, const std::vector<std::vector<int>>& weights) {
        max_weight_ = max_w;
        input_.resize(prices.size());
        for (int i = 0; i < input_.size(); ++i) {
            input_[i] = {prices[i], weights[i]};
        }
        cur_possible_sets_ = {{0, std::vector<int>(max_w.size())}};
    }
    
    bool CheckWeightSum(const WeightPrice& first, const WeightPrice& second) {
        for (int i = 0; i < second.weight.size(); ++i) {
            if (first.weight[i] + second.weight[i] > max_weight_[i]) {
                return false;
            }
        }
        return true;
    }

    WeightPrice Solve() {
        for (int i = 0; i < input_.size(); ++i) {
            std::set<WeightPrice> new_possible_sets;
            // std::cout << "AAAAAAAAAAAAAAAAA" << i << std::endl;
            for (const auto& state: cur_possible_sets_) {
                if (CheckWeightSum(state, input_[i])) {
                    new_possible_sets.insert(state + input_[i]);
                }
            }

            MergeLists(new_possible_sets);
        }
        return *--cur_possible_sets_.end();
    }

    void MergeLists(std::set<WeightPrice>& new_possible_sets) {
        std::set<WeightPrice> merged;
        auto first_ind = new_possible_sets.begin();
        auto second_ind = cur_possible_sets_.begin();
        std::vector<std::set<NDimensionalKnapsack::WeightPrice>::iterator> ind = {first_ind, second_ind};
        int picked = 0;
        if (*first_ind < *second_ind) {
            picked = 0;
        } else {
            picked = 1;
        }
        merged.insert(*ind[picked]);
        ++ind[picked];
        while (ind[0] != new_possible_sets.end() || ind[1] != cur_possible_sets_.end()) {
            if (ind[1]==cur_possible_sets_.end() || (ind[0] != new_possible_sets.end() && *ind[0] < *ind[1])) {
                picked = 0;
            } else {
                picked = 1;
            }
            if ((--merged.end())->price < ind[picked]->price) {
                if ((--merged.end())->weight == ind[picked]->weight || WeightPrice::VectorLess(ind[picked]->weight, (--merged.end())->weight)) {
                    merged.erase((--merged.end()));
                }
                merged.insert(*ind[picked]);
            }
            ++ind[picked];
        }
        cur_possible_sets_ = std::move(merged);
    }

private:
    std::vector<WeightPrice> input_;
    std::set<WeightPrice> cur_possible_sets_;
    std::vector<int> cur_picked_;
    std::vector<int> max_weight_;
};

// int main() {
//     NDimensionalKnapsack task;
//     task.Init({9}, {{6, {2}}, {5, {3}}, {8, {6}}, {9, {7}}, {6, {5}}, {7, {9}}, {3, {4}}});
//     auto ans = task.Solve();
//     std::cout << ans.price << " " << ans.weight[0] << std::endl;  
//     return 0;
// }
