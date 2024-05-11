#include <iostream>
#include <vector>
#include <set>


class OneDimensionalKnapsack {
public:
    struct WeightPrice {
        WeightPrice(int p, int w) : price(p), weight(w) {}

        int price = 0;
        int weight = 0;
        bool operator<(const  WeightPrice& other) const{
            return (this->price < other.price) || (this->price == other.price) && (this->weight < other.weight);
        }
        WeightPrice operator+(const  WeightPrice& other) const{
            return {this->price + other.price, this->weight + other.weight};
        }
    };

    OneDimensionalKnapsack() = default;

    void Init(int max_w, const std::vector<WeightPrice>& input) {
        max_weight_ = max_w;
        input_ = input;
        cur_possible_sets_ = {{0, 0}};
    }

    WeightPrice Solve() {
        for (int i = 0; i < input_.size(); ++i) {
            std::set<WeightPrice> new_possible_sets;
            std::cout << "AAAAAAAAAAAAAAAAA" << i << std::endl;
            for (const auto& state: cur_possible_sets_) {
                if (state.weight + input_[i].weight <= max_weight_) {
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
        std::vector<std::set<OneDimensionalKnapsack::WeightPrice>::iterator> ind = {first_ind, second_ind};
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
                if ((--merged.end())->weight >= ind[picked]->weight) {
                    merged.erase((--merged.end()));
                }
                merged.insert(*ind[picked]);
            }
            ++ind[picked];
        }
        for (auto a : new_possible_sets) {
            std::cout << a.price << " " << a.weight << std::endl;
        }
        std::cout << "NEW POSSIBLE" << std::endl;
        for (auto a : cur_possible_sets_) {
            std::cout << a.price << " " << a.weight << std::endl;
        }
        std::cout << "MERGED" << std::endl;
        for (auto a : merged) {
            std::cout << a.price << " " << a.weight << std::endl;
        }
        std::cout << "MERGED END" << std::endl;
        cur_possible_sets_ = std::move(merged);
    }

private:
    std::vector<WeightPrice> input_;
    std::set<WeightPrice> cur_possible_sets_;
    int max_weight_ = 0;
};

int main() {
    OneDimensionalKnapsack task;
    task.Init(9, {{6, 2}, {5, 3}, {8, 6}, {9, 7}, {6, 5}, {7, 9}, {3, 4}});
    auto ans = task.Solve();
    std::cout << ans.price << " " << ans.weight << std::endl;
    return 0;
}