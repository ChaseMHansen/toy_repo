import numpy as np
import matplotlib.pyplot as plt
import random


def generate_dataset(num_samples):
    """
    Generate a dataset from a Gaussian distribution with 0 mean and unit standard deviation
    """
    dataset = []
    for _ in range(num_samples):
        dataset.append(np.random.normal(loc=0.0, scale=1.0))
    return dataset


# P(A|B)=P(B|A)P(A)/P(B)


def p_x_given_a(x, a):
    N_a = np.sqrt(np.log(a) / (2 * np.pi))
    return N_a * a ** (-(x ** 2) / 2)


def pxa_many_x(x_list, a):
    """
    Find P({x}|a) from a list of x values.
    """
    result = 1
    for el in x_list:
        # The second term is to prevent the probability from going too low for numerics.
        # It is arbitrary to an extent. Ideally, I'd normalize, but I'm not sure exactly how.
        result *= p_x_given_a(el, a) * np.sqrt(4 * np.pi)
    return result


def metropolis(current_state, posterior, a_list=np.linspace(1, 4, 1000)):
    """
    Perform one step of the metropolis algorithm, does not move time forward.
    The generating function is just a random a value from the list of possible a values.
    """
    r = np.random.random()
    g = random.choice(range(len(a_list)))
    ratio = posterior[g] / posterior[current_state]
    if ratio >= 1:
        return g
    if ratio < r:
        return current_state
    if ratio > r:
        return g


def generate_posterior(num_samples, a_list=np.linspace(1, 4, 1000)):
    """
    Generate a posterior distribution from a random sample.
    """
    pax = []
    x_list = generate_dataset(num_samples)
    for a in a_list:
        # The lack of a second term is to keep P(a) uniform
        pax.append(pxa_many_x(x_list, a))

    return pax


def MCMC(num_iter, pax, a_list=np.linspace(1, 4, 1000)):
    """
    Run the Markov Chain Monte Carlo algorithm for num_iter steps on the posterior distribution.
    """
    posterior = pax
    current_state = random.choice(range(len(a_list)))
    chain = [a_list[current_state]]
    for _ in range(num_iter):
        link = metropolis(current_state, posterior, a_list)
        chain.append(a_list[link])
        current_state = link
    # Don't include the beginning of the chain to ensure a steady state.
    return chain[2000:]


pax = generate_posterior(1000)
chain = MCMC(5000, pax)
a_list = np.linspace(1, 4, 1000)
plt.subplot(3, 1, 1)
plt.plot(a_list, pax)
plt.subplot(3, 1, 2)
plt.plot(chain)
plt.subplot(3, 1, 3)
plt.hist(chain, bins=20)
plt.show()
