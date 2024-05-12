#include <vector>
#include <set>
#include <iostream>
#include <ctime>
#include <algorithm>
#include <cmath>
#include "ndim_branch_and_bound.cpp"


class OneDimensionalKnapsack {
public:

    OneDimensionalKnapsack() = default;

    WeightPrice Add(const WeightPrice& one, const WeightPrice& other) const{
        std::vector<int> tmp = one.weight;
        double tmp_price = one.price + other.price;
        for (int i = 1; i < one.weight.size(); ++i) {
            tmp[i] += other.weight[i];
            // std::cout << tmp_price << " " <<  (max_weights_[i] - tmp[i]) << std::endl;
            // tmp_price += lambdas_[i] * (max_weights_[i] - tmp[i]);
        }
        return {tmp_price, tmp};
    }

    void Init(const std::vector<int>& max_w, const std::vector<WeightPrice>& input, const std::vector<double>& lambdas) {
        max_weights_ = max_w;
        max_weight_ = max_w[0];
        input_ = input;
        lambdas_ = lambdas;
        for (int i = 0; i < input_.size(); ++i) {
            for (int j = 1; j < input_[i].weight.size(); ++j) {
                input_[i].price += lambdas_[j] * input_[i].weight[j];
            }
        }
        std::pair<WeightPrice, std::vector<int>> tmp = {{0, std::vector<int>(max_w.size(), 0)}, std::vector<int>()};
        cur_possible_sets_ = {tmp};
    }

    std::pair<WeightPrice, std::vector<int>> Solve() {
        for (int i = 0; i < input_.size(); ++i) {
            std::set<std::pair<WeightPrice, std::vector<int>>> new_possible_sets;
            for (const auto& state: cur_possible_sets_) {
                if (state.first.weight[0] + input_[i].weight[0] <= max_weight_) {
                    auto tmp = state.second;
                    tmp.push_back(i);
                    new_possible_sets.insert({Add(state.first, input_[i]), tmp});
                }
            }
            MergeLists(new_possible_sets);
        }
        return *--cur_possible_sets_.end();
    }

    void MergeLists(std::set<std::pair<WeightPrice, std::vector<int>>>& new_possible_sets) {
        std::set<std::pair<WeightPrice, std::vector<int>>> merged;
        auto first_ind = new_possible_sets.begin();
        auto second_ind = cur_possible_sets_.begin();
        std::vector<std::set<std::pair<WeightPrice, std::vector<int>>>::iterator> ind = {first_ind, second_ind};
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
            if ((--merged.end())->first.price < ind[picked]->first.price) {
                if ((--merged.end())->first.weight[0] >= ind[picked]->first.weight[0]) {
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
    std::set<std::pair<WeightPrice, std::vector<int>>> cur_possible_sets_;
    std::vector<int> max_weights_;
    int max_weight_ = 0;
    std::vector<double> lambdas_;
};


class NDimensionalLagrange {
public:
    void Init(const std::vector<int>& max_w, const std::vector<int>& prices, const std::vector<std::vector<int>>& weights) {
        max_weight_ = max_w;
        // input_ = input;
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
        ++n;
        // if (count % 1000 == 0) {
        //     std::cout << count << std::endl;
        // }
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
        std::cout << SolveLagrange(node - 1) + cur_price_weight_.price << " " << best_price_ << std::endl;
        if (best_price_ == cur_price_weight_.price || SolveLagrange(node - 1) + cur_price_weight_.price > best_price_) {
            // std::cout << "NO BLOCK" << std::endl;
            Dfs(node - 1, true);
            Dfs(node - 1, false);
        }
        // } else {
        //     std::cout << "blocked" << std::endl;
        // }
        if (pick) {
            cur_price_weight_ = cur_price_weight_ - input_[node];
            cur_picked_.pop_back();
        }
    }

    double LagrangianObj(const std::vector<int>& x, double price, const std::vector<double>& lambdas) {
        double tmp = price;
        std::vector<int> subgrad = max_weight_;
        for (auto elem : x) {
            for (int i = 1; i < max_weight_.size(); ++i) {
                subgrad[i] = input_[elem].weight[i];
            }
        }
        for (int i = 1; i < lambdas.size(); ++i) {
            tmp += lambdas[i] * subgrad[i];
        }
        return tmp;
    }

    double SolveLagrange(int nodes_left) {
        OneDimensionalKnapsack relax;
        auto tmp = max_weight_;
        for (int i = 0; i < tmp.size(); ++i) {
            tmp[i] -= cur_price_weight_.weight[i];
        }
        auto lambdas = std::vector<double>(max_weight_.size(), 1);
        std::pair<WeightPrice, std::vector<int>> ans;
        double last_ans = 0;
        double ret = 0;
        std::cout << "NEW";
        for (int i = 0; i < 10; ++i) {
            relax.Init(tmp, std::vector<WeightPrice>(input_.begin(), input_.begin() + nodes_left + 1), lambdas);
            ans = relax.Solve();
            if (std::abs(ans.first.price - last_ans) < 0.1) {
                break;
            }
            ret = LagrangianObj(ans.second, ans.first.price, lambdas);
            Descent(lambdas, ans.second);
        }
        return ret;
    }

    void Descent(std::vector<double>& lambdas, const std::vector<int>& indices) {
        std::vector<int> subgrad = max_weight_;
        subgrad[0] = 0;
        for (auto elem : indices) {
            for (int i = 1; i < max_weight_.size(); ++i) {
                subgrad[i] -= input_[elem].weight[i];
            }
        }
        double norm = 0;
        for (auto a : subgrad) {
            norm += a * a;
        }
        for (int i = 1; i < lambdas.size(); ++i) {
            lambdas[i] -= 1.0/std::sqrt(norm) * subgrad[i];
            lambdas[i] = std::max<double>(0, lambdas[i]);
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


std::vector<std::vector<int>> GenerateInputWeights(int ndim, int input_size, int mod) {
    std::vector<std::vector<int>> gen_input(input_size, std::vector<int>(ndim));
    for (int i = 0; i < input_size; ++i) {
        for (int j = 0; j < ndim; ++j) {
            gen_input[i][j] = std::rand() % mod + 1;
        }
    }
    return gen_input;
}

std::vector<int> GenerateInputPrices(int ndim, int input_size, int mod) {
    std::vector<int> gen_input(input_size);
    for (int i = 0; i < input_size; ++i) {
        gen_input[i] = std::rand() % mod + 1;
    }
    return gen_input;
}

std::vector<int> GenerateWeight(int ndim, int mod) {
    std::vector<int> gen_max_w(ndim);
    std::generate(gen_max_w.begin(), gen_max_w.end(), std::rand);
    for (int i = 0; i < ndim; ++i) {
        gen_max_w[i] = gen_max_w[i] % mod + 1+ mod;
    }
    return gen_max_w;
}

int main() {
    int iterations = 1;
    int ndim = 20;
    int n_items = 10;
    std::vector<double> timers_lag(iterations);
    // std::vector<double> timers_dp(iterations);

    std::vector<double> timers_bb(iterations);
    std::vector<int> n_bb(iterations);
    std::vector<int> n_lag(iterations);
    for (int i = 0; i < iterations; ++i) {
        std::cout <<"Iteration " << i << std::endl;
        NDimensionalLagrange task;
        NDimensionalBranch task3;
        std::srand(unsigned(std::time(nullptr)+i));
        auto gen_w = GenerateWeight(ndim, 50);
        auto gen_input_weight = GenerateInputWeights(ndim, n_items, 20);
        auto gen_input_price = GenerateInputPrices(ndim, n_items, 20);
        task.Init(gen_w, gen_input_price, gen_input_weight);
        task3.Init(gen_w, gen_input_price, gen_input_weight);
        // std::cout<< "INIT DONE" << std::endl;
        clock_t start = clock();
        auto ans1 = task.Solve();
        n_lag[i] = task.n;
        clock_t end = clock();
        double seconds = (double)(end - start) / CLOCKS_PER_SEC;
        timers_lag[i] = seconds;
        start = clock();
        auto ans3 = task3.Solve();
        // std::cout << ans3 << std::endl;
        end = clock();
        seconds = (double)(end - start) / CLOCKS_PER_SEC;
        // std::cout << seconds << std::endl;
        timers_bb[i] = seconds;
        n_bb[i] = task3.n;
    }
    for (int i = 0; i < iterations; ++i) {
        std::cout << timers_lag[i] << " ";
    }
    std::cout << std::endl;
    for (int i = 0; i < iterations; ++i) {
        std::cout << timers_bb[i] << " ";
    }
    std::cout << std::endl;
        for (int i = 0; i < iterations; ++i) {
        std::cout << n_lag[i] << " ";
    }
    std::cout << std::endl;
    for (int i = 0; i < iterations; ++i) {
        std::cout << n_bb[i] << " ";
    }
    return 0;
}
