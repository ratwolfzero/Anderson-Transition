import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
import qutip as qt
from dataclasses import dataclass
import logging
from scipy.sparse import diags
import time

# -----------------------
# Logging
# -----------------------
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# -----------------------
# Simulation parameters
# -----------------------


@dataclass
class SimulationParams:
    N: int = None
    t: float = None
    disorder_strength: float = None
    dt: float = 0.5
    total_time: float = None
    sigma: float = None
    seed: int = None

    def __post_init__(self):
        if self.N <= 0 or self.dt <= 0 or self.total_time <= 0:
            raise ValueError("N, dt, and total_time must be positive")
        if self.sigma <= 0:
            raise ValueError("Gaussian width (sigma) must be positive")

    @property
    def frames(self) -> int:
        return int(self.total_time / self.dt)


def get_default_params() -> SimulationParams:
    return SimulationParams(
        N=100,
        t=1.0,
        disorder_strength=2.5,
        dt=0.3,
        total_time=25,
        sigma=4.0,
        seed=42
    )

# -----------------------
# Initial state
# -----------------------


def create_initial_state(sites: np.ndarray, params: SimulationParams, state_type: str = 'gaussian') -> qt.Qobj:
    if state_type == 'gaussian':
        return qt.Qobj(np.exp(-(sites - params.N//2)**2 / (2 * params.sigma**2))).unit()
    elif state_type == 'localized':
        psi = np.zeros(params.N)
        psi[params.N//2] = 1.0
        return qt.Qobj(psi).unit()
    raise ValueError(f"Unknown state type: {state_type}")

# -----------------------
# Hamiltonian
# -----------------------


def build_hamiltonian(params: SimulationParams, disorder_strength: float = None) -> qt.Qobj:
    N, t = params.N, params.t
    W = disorder_strength if disorder_strength is not None else params.disorder_strength
    hop = diags([np.ones(N-1), np.ones(N-1)], offsets=[1, -1], shape=(N, N))
    disorder = W * (np.random.rand(N) - 0.5)
    onsite = diags(disorder, shape=(N, N))
    return qt.Qobj(-t * hop + onsite, dims=[[N], [N]])

# -----------------------
# Time evolution
# -----------------------


def evolve_probs(H: qt.Qobj, psi0: qt.Qobj, times: np.ndarray) -> np.ndarray:
    result = qt.mesolve(H, psi0, times, options={'progress_bar': 'enhanced'})
    return np.abs(np.array([state.full().flatten() for state in result.states]))**2

# -----------------------
# Participation ratio
# -----------------------


def participation_ratio(probs: np.ndarray) -> np.ndarray:
    return 1.0 / np.sum(probs**2, axis=1)

# -----------------------
# Optimized Animation with dynamic PR and legends
# -----------------------


def create_animation(fig, ax_clean, ax_disordered, ax_pr_clean, ax_pr_disordered,
                     sites: np.ndarray, probs_clean: np.ndarray, probs_disordered: np.ndarray,
                     times: np.ndarray, pr_clean: np.ndarray, pr_disordered: np.ndarray) -> FuncAnimation:

    # Clean system: filled area
    fill_clean = ax_clean.fill_between(
        sites, probs_clean[0], color='blue', alpha=0.3)
    # Disordered system: line only
    line_disordered, = ax_disordered.plot(
        sites, probs_disordered[0], color='red', alpha=0.8)

    time_text = ax_disordered.text(0.95, 0.95, '', transform=ax_disordered.transAxes,
                                   ha='right', va='top', fontsize=10)

    # Participation ratio lines (initially empty)
    pr_line_clean, = ax_pr_clean.plot([], [], color='blue')  # Removed label
    pr_line_disordered, = ax_pr_disordered.plot(
        [], [], color='red')  # Removed label

    # Removed initial legend call

    # Set axes limits
    ax_pr_clean.set_xlim(times[0], times[-1])
    ax_pr_clean.set_ylim(0, np.max(pr_clean)*1.1)
    ax_pr_disordered.set_xlim(times[0], times[-1])
    ax_pr_disordered.set_ylim(0, np.max(pr_disordered)*1.1)

    def init():
        nonlocal fill_clean
        fill_clean.remove()
        fill_clean = ax_clean.fill_between(
            sites, np.zeros_like(sites), color='blue', alpha=0.3)
        line_disordered.set_ydata(np.zeros_like(sites))
        time_text.set_text('')
        pr_line_clean.set_data([], [])
        pr_line_disordered.set_data([], [])
        return fill_clean, line_disordered, time_text, pr_line_clean, pr_line_disordered

    def animate(i):
        nonlocal fill_clean
        # Update probability densities
        fill_clean.remove()
        fill_clean = ax_clean.fill_between(
            sites, probs_clean[i], color='blue', alpha=0.3)
        line_disordered.set_ydata(probs_disordered[i])
        time_text.set_text(f"t = {times[i]:.1f}")

        # Update PR lines
        pr_line_clean.set_data(times[:i+1], pr_clean[:i+1])
        pr_line_disordered.set_data(times[:i+1], pr_disordered[:i+1])

        # Removed legend calls

        return fill_clean, line_disordered, time_text, pr_line_clean, pr_line_disordered

    return FuncAnimation(fig, animate, frames=len(times), init_func=init, blit=True)

# -----------------------
# Visualization setup
# -----------------------


def setup_visualization(sites, probs_clean, probs_disordered, times,
                        pr_clean, pr_disordered, params):
    plt.style.use('ggplot')
    fig, axs = plt.subplots(2, 2, figsize=(
        12, 8), gridspec_kw={'height_ratios': [2, 1]})

    max_clean = np.max(probs_clean)
    max_disordered = np.max(probs_disordered)

    # Top row: probability density
    axs[0, 0].set_title('Clean System (Delocalized)')
    axs[0, 1].set_title(
        f'Disordered System (W={params.disorder_strength}, Localized)')
    axs[0, 0].set_xlabel('Lattice Site')
    axs[0, 1].set_xlabel('Lattice Site')
    axs[0, 0].set_ylabel('Probability Density |ψ(x,t)|²')
    for i, ax in enumerate(axs[0]):
        ax.set_xlim(0, params.N-1)
        ax.set_ylim(0, max_clean*1.1 if i == 0 else max_disordered*1.1)
        ax.axvline(params.N//2, color='gray', ls='--', lw=1)

    # Bottom row: participation ratio axes
    axs[1, 0].set_title("Participation Ratio (Clean)")
    axs[1, 1].set_title("Participation Ratio (Disordered)")
    axs[1, 0].set_xlabel("Time")
    axs[1, 1].set_xlabel("Time")
    axs[1, 0].set_ylabel("PR")
    axs[1, 0].set_xlim(0, params.total_time)
    axs[1, 1].set_xlim(0, params.total_time)
    axs[1, 0].set_ylim(0, np.max(pr_clean)*1.1)
    axs[1, 1].set_ylim(0, np.max(pr_disordered)*1.1)
    axs[1, 0].grid(True)
    axs[1, 1].grid(True)

    # Animation with dynamic PR and legends
    ani = create_animation(fig, axs[0, 0], axs[0, 1], axs[1, 0], axs[1, 1],
                           sites, probs_clean, probs_disordered, times, pr_clean, pr_disordered)
    return fig, axs, ani

# -----------------------
# Main simulation
# -----------------------


def main(params: SimulationParams, save_animation=False, save_plot=False):
    logger.info("Starting simulation with N=%d, W=%.2f",
                params.N, params.disorder_strength)
    np.random.seed(params.seed)
    start_init = time.time()

    sites = np.arange(params.N)
    times = np.linspace(0, params.total_time, params.frames)
    psi0 = create_initial_state(sites, params)

    H_clean = build_hamiltonian(params, disorder_strength=0.0)
    H_disordered = build_hamiltonian(params)
    probs_clean = evolve_probs(H_clean, psi0, times)
    probs_disordered = evolve_probs(H_disordered, psi0, times)

    pr_clean = participation_ratio(probs_clean)
    pr_disordered = participation_ratio(probs_disordered)

    logger.info("Initialization done in %.4f s", time.time() - start_init)

    start_viz = time.time()
    fig, axs, ani = setup_visualization(
        sites, probs_clean, probs_disordered, times, pr_clean, pr_disordered, params)
    plt.tight_layout()

    if save_animation:
        ani.save("anderson_localization.gif",
                 writer=PillowWriter(fps=10), dpi=150)
    if save_plot:
        plt.savefig("participation_ratio.png", dpi=150)

    plt.show()
    logger.info("Visualization done in %.4f s", time.time() - start_viz)


# -----------------------
# Run simulation
# -----------------------
if __name__ == "__main__":
    params = get_default_params()
    main(params, save_animation=False, save_plot=False)
