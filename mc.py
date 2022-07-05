from collections import namedtuple
from math import e, exp  
from random import random


GenState = namedtuple(
    "GenState",
    [
        "uninfected_bacteria",
        "phages",
        "infected_bacteria",
        "carry_over",
        "gen_count",
    ],
)
MCParameters = namedtuple(
    "MCParameters",
    [
        "growth_rate",
        "carrying_capacity",
        "burst_size",
        "latent_period",
        "secondary_decay_proportionality_constant",
    ],
)


def calculate_next_gen(CurrentState, SimulationParameters):
    (
        uninfected_bacteria,
        phages,
        infected_bacteria,
        carry_over,
        _,
    ) = CurrentState

    increased_uninfected_bacteria = (
        SimulationParameters.growth_rate
        * uninfected_bacteria
        * (1 - uninfected_bacteria / SimulationParameters.carrying_capacity)
    )
    dN = increased_uninfected_bacteria - uninfected_bacteria
    duplication_rate = dN / uninfected_bacteria

    moi = phages / uninfected_bacteria
    adsorption_rate = moi * exp(-moi)

    carry_over_index = CurrentState.gen_count + SimulationParameters.latent_period
    carry_over_size = len(carry_over)

    carry_over[carry_over_index % carry_over_size] = {
        "phage_count": 0,
        "infected_bacteria_count": 0,
        "uninfected_bacteria_count": 0,
    }

    for _ in range(uninfected_bacteria):
        random_1 = random()
        random_2 = random()
        if random_1 >= adsorption_rate:
            if random_2 >= duplication_rate:
                continue
            elif random_2 < duplication_rate:
                uninfected_bacteria += 1
        elif random_1 < adsorption_rate:
            uninfected_bacteria -= 1
            infected_bacteria += 1
            phages -= 1

            carry_over[carry_over_index % carry_over_size][
                "phage_count"
            ] += SimulationParameters.burst_size
            carry_over[carry_over_index % carry_over_size][
                "infected_bacteria_count"
            ] -= 1

            secondary_killing = int(
                uninfected_bacteria
                * exp(
                    -CurrentState.gen_count
                    / SimulationParameters.latent_period
                    * SimulationParameters.secondary_decay_proportionality_constant
                )
            )
            carry_over[carry_over_index % carry_over_size][
                "uninfected_bacteria_count"
            ] -= secondary_killing

    carry_over[carry_over_index % carry_over_size]["phage_count"] //= 20

    NextGen = GenState(
        uninfected_bacteria,
        phages,
        infected_bacteria,
        carry_over,
        CurrentState.gen_count,
    )
    return NextGen


def apply_carry_over(NextGen):
    phages = NextGen.phages
    infected_bacteria = NextGen.infected_bacteria
    uninfected_bacteria = NextGen.uninfected_bacteria
    gen_count = NextGen.gen_count

    carry_over_size = len(NextGen.carry_over)

    phages += NextGen.carry_over[gen_count % carry_over_size]["phage_count"]
    infected_bacteria += NextGen.carry_over[gen_count % carry_over_size][
        "infected_bacteria_count"
    ]

    secondary_killing_rate = abs(
        NextGen.carry_over[gen_count % carry_over_size]["uninfected_bacteria_count"]
        / uninfected_bacteria
    )
    for _ in range(uninfected_bacteria):
        r = random()
        if r < secondary_killing_rate:
            uninfected_bacteria -= 1

    gen_count += 1
    UpdatedNextGen = GenState(
        NextGen.uninfected_bacteria,
        phages,
        infected_bacteria,
        NextGen.carry_over,
        gen_count,
    )
    return UpdatedNextGen


def monte_carlo(n, InitialState, SimulationParameters):
    for _ in range(n):
        InitialState = apply_carry_over(
            calculate_next_gen(InitialState, SimulationParameters)
        )
    return InitialState

def get_mc_input():
    def _get_valid_input(input_prompt, wrong_input_prompt, is_valid):
        answer = input(input_prompt)
        while not is_valid(answer):
            print(wrong_input_prompt)
            answer = input(input_prompt)
        return answer

    UNINFECTED_BACTERIA_LOWER_BOUND = 1000
    UNINFECTED_BACTERIA_UPPER_BOUND = 10000

    IsValid = namedtuple("IsValid", ["uninfected_bacteria", "moi", "growth_rate", "burst_size", "latent_period", "monte_carlo_times"])
    MCIsValid = IsValid(
    lambda x: int(x) in range(UNINFECTED_BACTERIA_LOWER_BOUND, UNINFECTED_BACTERIA_UPPER_BOUND + 1),
    lambda x: float(x) in [.1, 1],
    lambda x: float(x) >= .02 and float(x) <= 1,
    lambda x: int(x) in range(100, 301),
    lambda x: int(x) in range(2, 11),
    lambda x: int(x) in range(1, 10001))

    Prompt = namedtuple("Prompt", ["uninfected_bacteria", "moi", "growth_rate", "carrying_capacity", "burst_size", "latent_period", "monte_carlo_times"])
    MCPrompt = Prompt(
        ("Uninfected bacteria (U) [1000, 10000]: \n", "Must be in [1000, 10000]"),
        ("Multiplicity of infection (MOI) (P/U) [0.1 | 1]: \n", "Must be either 0.1 or 1"),
        ("Growth rate (μ) [0.2, 1]: \n", "Must be in [0.2, 1]"),
        ("Carrying capacity (Κ) [2*U, 10*U]: \n", "Must be in [2*U, 10*U]"),
        ("Burst size (bS) [100, 300]: \n", "Must be in [100, 300]"),
        ("Latent period (l) [2, 10]: \n", "Must be in [2, 10*e]"),
        ("Monte Carlo simulation times [1, 10000]: \n", "Must be in [1, 10000]")
    )

    uninfected_bacteria = int(_get_valid_input(MCPrompt.uninfected_bacteria[0], MCPrompt.uninfected_bacteria[1], MCIsValid.uninfected_bacteria))
    moi = float(_get_valid_input(MCPrompt.moi[0], MCPrompt.moi[1], MCIsValid.moi))
    growth_rate = float(_get_valid_input(MCPrompt.growth_rate[0], MCPrompt.growth_rate[1], MCIsValid.growth_rate))

    is_valid_carrying_capacity = lambda x: int(x) >= 2*uninfected_bacteria and int(x) <= 10*uninfected_bacteria
    carrying_capacity = int(_get_valid_input(MCPrompt.carrying_capacity[0], MCPrompt.carrying_capacity[1], is_valid_carrying_capacity))

    burst_size = int(_get_valid_input(MCPrompt.burst_size[0], MCPrompt.burst_size[1], MCIsValid.burst_size))
    latent_period = int(_get_valid_input(MCPrompt.latent_period[0], MCPrompt.latent_period[1], MCIsValid.latent_period))
    monte_carlo_times = int(_get_valid_input(MCPrompt.monte_carlo_times[0], MCPrompt.monte_carlo_times[1], MCIsValid.monte_carlo_times))

    return uninfected_bacteria, moi, growth_rate, carrying_capacity, burst_size, latent_period, monte_carlo_times

def setup(user_input):
    uninfected_bacteria, moi, growth_rate, carrying_capacity, burst_size, latent_period, monte_carlo_times = user_input

    phages = int(uninfected_bacteria * moi)
    infected_bacteria, gen_count = 0, 0

    carry_over = [
        {"phage_count": 0, "infected_bacteria_count": 0, "uninfected_bacteria_count": 0}
        for _ in range(latent_period + 1)
    ]
    secondary_decay_proporitionality_constant = 2 if moi == .1 else e/2

    InitialState = GenState(uninfected_bacteria, phages, infected_bacteria, carry_over, gen_count)
    SimulationParameters = MCParameters(growth_rate, carrying_capacity, burst_size, latent_period, secondary_decay_proporitionality_constant)

    return monte_carlo_times, InitialState, SimulationParameters


def main():
    n, InitialState, SimulationParameters = setup(get_mc_input())
    FinalState = monte_carlo(n, InitialState, SimulationParameters)

    print(f"Uninfected bacteria count: {FinalState.uninfected_bacteria}")
    print(f"Phages count: {FinalState.phages}")
    print(f"Infected bacteria: {FinalState.infected_bacteria}")

main()