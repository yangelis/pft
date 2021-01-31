#include "../pft.hpp"
#include <random>
#include <unistd.h>

enum class State : i32 { Dead = 0, Alive = 1 };

using Cell = pft::Matrix<State>;

struct Board {
  Cell cells;

  auto wrap(i32 row, i32 col) const -> State {
    return cells(pft::mod(row, (i32)cells.rows),
                 pft::mod(col, (i32)cells.cols));
  }

  auto count_neighbors(i32 y, i32 x) const -> i32 {
    i32 acc = 0;
    for (i32 dy = -1; dy <= 1; ++dy) {
      for (i32 dx = -1; dx <= 1; ++dx) {
        if (dx != 0 || dy != 0) {
          if (wrap(y + dy, x + dx) == State::Alive) {
            acc += 1;
          }
        }
      }
    }
    return acc;
  }

  void render() const {
    for (std::size_t i = 0; i < cells.rows; ++i) {
      for (std::size_t j = 0; j < cells.cols; ++j) {
        if (cells(i, j) == State::Alive) {
          pft::print(stdout, '#');
        } else {
          pft::print(stdout, '.');
        }
      }
      pft::print(stdout, '\n');
    }
  }
};

void next_gen(const Board& prev, Board& next) {
  for (size_t i = 0; i < prev.cells.rows; ++i) {
    for (size_t j = 0; j < prev.cells.cols; ++j) {
      const auto nbors = prev.count_neighbors(i, j);
      if (prev.cells(i, j) == State::Alive) {
        next.cells(i, j) =
            nbors == 3 || nbors == 2 ? State::Alive : State::Dead;
      } else {
        next.cells(i, j) = nbors == 3 ? State::Alive : State::Dead;
      }
    }
  }
}

int main(int argc, char* argv[]) {
  Board boards[2];
  size_t Width  = 10;
  size_t Height = 10;

  if (argc == 1) {
    boards[0].cells = Cell{Height, Width};
    boards[1].cells = Cell{Height, Width};
    // set up the glider
    boards[0].cells(1, 2) = State::Alive;
    boards[0].cells(2, 3) = State::Alive;
    boards[0].cells(3, 1) = State::Alive;
    boards[0].cells(3, 2) = State::Alive;
    boards[0].cells(3, 3) = State::Alive;
  } else if (argc == 3) {
    Width           = pft::to_int(argv[1]);
    Height          = pft::to_int(argv[2]);
    boards[0].cells = Cell{Height, Width};
    boards[1].cells = Cell{Height, Width};
    // randomize the state
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, 1);
    for (size_t i = 0; i < Height; ++i) {
      for (size_t j = 0; j < Width; ++j) {
        boards[0].cells(i, j) = State(dis(gen));
      }
    }
  } else {
    pft::panic("Bad arguments");
  }

  i32 ii = 0;

  for (;;) {
    boards[ii].render();
    i32 jj = 1 - ii;
    next_gen(boards[ii], boards[jj]);
    ii = jj;
    pft::print(stdout, "\033[", Height, "A");
    pft::print(stdout, "\033[", Width, "D");
    usleep(700000);
  }
  return 0;
}
