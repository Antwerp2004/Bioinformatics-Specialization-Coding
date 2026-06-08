import sys
import math
from datatypes import Cell, Board


def initialize_board(num_rows: int, num_cols: int) -> Board:
    return [[(0.0, 0.0) for _ in range(num_cols)] for _ in range(num_rows)]


def sum_cells(*cells: Cell) -> list[float]:
    """
    Sum corresponding elements of multiple cells.

    Args:
        *cells: An arbitrary number of Cell values, where each Cell is
            a tuple of two floats.
    Returns:
        A single Cell representing the element-wise sum of all input cells.
    """
    sum_cell_prey, sum_cell_predator = 0, 0
    for prey, predator in cells:
        sum_cell_prey += prey
        sum_cell_predator += predator
    return sum_cell_prey, sum_cell_predator


def change_due_to_reactions(
    current_cell: Cell,
    feed_rate: float,
    kill_rate: float
) -> Cell:
    """
    Compute the change in a cell due to Gray-Scott reactions.

    Args:
        current_cell: The current cell.
        feed_rate: The feed reaction rate.
        kill_rate: The kill reaction rate.
    Returns:
        A Cell representing the change in current_cell due to reactions.
    """
    prey, predator = current_cell
    prey, predator = feed_rate * (1 - prey) - prey * predator**2, - kill_rate * predator + prey * predator**2
    return prey, predator


def change_due_to_diffusion(
    current_board: Board,
    row: int,
    col: int,
    prey_diffusion_rate: float,
    predator_diffusion_rate: float,
    kernel: list[list[float]]
) -> Cell:
    """
    Compute the change in a cell due to diffusion.

    Args:
        current_board: A 2D list of Cells, where each Cell is a tuple of two floats.
        row: Row index of the cell.
        col: Column index of the cell.
        prey_diffusion_rate: Diffusion rate for the prey component.
        predator_diffusion_rate: Diffusion rate for the predator component.
        kernel: A 3x3 diffusion kernel.
    Returns:
        A Cell representing the change in the cell at (row, col) due to diffusion.
    """
    
    rows = len(current_board)
    cols = len(current_board[0]) if rows > 0 else 0
    prey, predator = 0, 0
    
    for i in range(-1, 2):
        for j in range(-1, 2):
            if (0 <= row + i <= rows - 1) and (0 <= col + j <= cols - 1):
                prey_change_rate = prey_diffusion_rate * kernel[i + 1][j + 1]
                predator_change_rate = predator_diffusion_rate * kernel[i + 1][j + 1]
                    
                prey += prey_change_rate * current_board[row + i][col + j][0] 
                predator += predator_change_rate * current_board[row + i][col + j][1]
     
    return prey, predator


def update_cell(current_board: Board, row: int, col: int,
                feed_rate: float, kill_rate: float,
                prey_diffusion_rate: float, predator_diffusion_rate: float, kernel: list[list[float]]):
    current_cell = current_board[row][col]
    diffusion_values = change_due_to_diffusion(current_board, row, col, prey_diffusion_rate, predator_diffusion_rate, kernel)
    reaction_values = change_due_to_reactions(current_cell, feed_rate, kill_rate)
    return sum_cells(current_cell, diffusion_values, reaction_values)


def update_board(current_board: Board, feed_rate: float, kill_rate: float,
                 prey_diffusion_rate: float, predator_diffusion_rate: float, kernel: list[list[float]]) -> Board:
    """
    Update a Gray-Scott reaction-diffusion board by one time step.

    Args:
        current_board: A 2D list of Cells, where each Cell is a list of two floats.
        feed_rate: Feed reaction rate.
        kill_rate: Kill reaction rate.
        prey_diffusion_rate: Diffusion rate for the prey component.
        predator_diffusion_rate: Diffusion rate for the predator component.
        kernel: A 3x3 diffusion kernel.
    Returns:
        A new Board representing the next time step after applying the
        Gray-Scott reaction-diffusion update rules.
    """
    num_rows = len(current_board)
    num_cols = len(current_board[0]) if num_rows > 0 else 0
    new_board = [[0] * num_cols for _ in range(num_rows)]
    for row in range(num_rows):
        for col in range(num_cols):
            new_board[row][col] = update_cell(current_board, row, col, feed_rate, kill_rate, prey_diffusion_rate, predator_diffusion_rate, kernel)
    return new_board


def simulate_gray_scott(
    initial_board: Board,
    num_gens: int,
    feed_rate: float,
    kill_rate: float,
    prey_diffusion_rate: float,
    predator_diffusion_rate: float,
    kernel: list[list[float]]
) -> Board:
    """
    Simulate the Gray-Scott reaction-diffusion system.

    Args:
        initial_board: A 2D list of Cells, where each Cell is a list of two floats.
        num_gens: Number of generations (time steps) to simulate.
        feed_rate: Feed reaction rate.
        kill_rate: Kill reaction rate.
        prey_diffusion_rate: Diffusion rate for the prey component.
        predator_diffusion_rate: Diffusion rate for the predator component.
        kernel: A 3x3 diffusion kernel.
    Returns:
        A list of Boards of length num_gens + 1, where the first board is
        initial_board and each subsequent board is the next time step.
    """
    boards = [0] * (num_gens + 1)
    boards[0] = initial_board
    for i in range(1, num_gens + 1):
        boards[i] = update_board(boards[i-1], feed_rate, kill_rate, prey_diffusion_rate, predator_diffusion_rate, kernel)
    return boards