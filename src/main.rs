use std::time::Duration;
use std::ops::{Index,IndexMut,Mul,Add,AddAssign};
use rand::Rng;

use ndarray::{array,Array1,Array2};

use sdl2::event::Event;
use sdl2::keyboard::Keycode;
use sdl2::pixels::Color;
use sdl2::render::Canvas;
use sdl2::video::Window;
use sdl2::gfx::primitives::DrawRenderer;

const NDF: usize = 10;
const DIM: usize = 2*NDF;

fn sqr(x: f64) -> f64 {
    x*x
}

#[derive(Clone, Copy)]
struct State {
    f: [f64; DIM],
}

impl Index<usize> for State {
    type Output = f64;
    fn index(&self, idx: usize) -> &f64 {
        &self.f[idx]
    }
}
impl IndexMut<usize> for State {
    fn index_mut(&mut self, idx: usize) -> &mut f64 {
        &mut self.f[idx]
    }
}
impl Mul<f64> for State {
    type Output = Self;
    fn mul(self, rhs: f64) -> Self {
        let mut result = Self::zero();
        for i in 0..DIM {
            result[i] = self[i] * rhs;
        }
        result
    }
}
impl Mul<State> for f64 {
    type Output = State;
    fn mul(self, rhs: State) -> State {
        rhs*self
    }
}
impl Add for State {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        let mut result = Self::zero();
        for i in 0..DIM { result[i] = self[i] + rhs[i] };
        result
    }
}
impl AddAssign for State {
    fn add_assign(&mut self, rhs: Self) {
        for i in 0..DIM { self[i] += rhs[i] };
    }
}

impl State {
    fn zero() -> State {
        State{ f: [0.; DIM] }
    }
    fn from_halves(f1: &[f64], f2: [f64; NDF]) -> State {
        let mut f = [0.; DIM];
        for i in 0..NDF {
            f[i] = f1[i];
            f[NDF+i] = f2[i];
        }
        State {f}
    }
    fn from_rand() -> State {
        let mut rng = rand::thread_rng();
        let mut f = [0.; DIM];
        for i in 0..NDF {
            f[i] = 8.*rng.gen::<f64>()-0.5;
            f[NDF+i] = 4.*rng.gen::<f64>()-0.5;
        }
        State { f }
    }

    fn deriv(p: State, dir: i32) -> State {
        // s_{jk} = min(j+1, k+1)*sin(phi_j-phi_k)
        // c_{jk} = min(j+1, k+1)*cos(phi_j-phi_k)
        let mut s = [[0.; NDF]; NDF];
        let mut c = [[0.; NDF]; NDF];
        for i in 0..NDF {
            for j in 0..i {
                s[i][j] = (j+1) as f64*(p[i]-p[j]).sin();
                s[j][i] = -s[i][j];
                c[i][j] = (j+1) as f64*(p[i]-p[j]).cos();
                c[j][i] = c[i][j];
            }
        }
        // RHS vector
        let mut y = [0.; NDF];
        for i in 0..NDF {
            y[i] = 0.5 * i as f64 * p[i].sin();
            for j in 0..NDF {
                y[i] += s[i][j] * sqr(p[NDF+j]);
            }
        }
        // TODO: Invert c*x=y using CG

        State::from_halves(
            &p.f[NDF..DIM],
            [0.; NDF],
        )
    }

    fn energy(&self) -> f64 {
        0.
    }

    fn step(&mut self, dir: i32) {
        let h = 0.001;
        let k1 = Self::deriv(*self, dir);
        let k2 = Self::deriv(*self + 0.5*h*k1, dir);
        let k3 = Self::deriv(*self + 0.5*h*k2, dir);
        let k4 = Self::deriv(*self + h*k3, dir);
        *self += h/6.0 * (k1+2.*k2+2.*k3+k4);
    }

    fn draw(&self, canvas: &mut Canvas<Window>) -> Result<(), String> {
        println!("{}", self.energy());
        let mut pos = (512., 384.);
        let mut nextpos: (f64, f64);
        for i in (0..NDF).rev() {
            let (s, c) = self[i].sin_cos();
            let nextpos = (
                pos.0 + 50. * s,
                pos.1 + 50. * c,
            );
            canvas.thick_line(
                pos.0 as i16, pos.1 as i16,
                nextpos.0 as i16, nextpos.1 as i16,
                10, Color::RGB(0,0,0),
            )?;
            canvas.filled_circle(
                pos.0 as i16, pos.1 as i16, 5, Color::RGB(0,0,0)
            )?;
            pos = nextpos;
        }
        Ok(())
    }
}

fn main() -> Result<(), String> {
    let sdl_context = sdl2::init()?;
    let video = sdl_context.video()?;
    let width = 1024;
    let height = 800;
    let mut canvas = video
        .window("Pendulum", width, height)
        .position_centered()
        .build()
        .map_err(|e| e.to_string())?
        .into_canvas()
        .accelerated()
        .build()
        .map_err(|e| e.to_string())?;
    let mut event_pump = sdl_context.event_pump()?;

    let mut state = State::from_rand();
    let mut dir = 0;


    'running: loop {
        canvas.set_draw_color(Color::RGB(155, 155, 155));
        canvas.clear();
        for event in event_pump.poll_iter() {
            match event {
                Event::Quit {..} |
                Event::KeyDown { keycode: Some(Keycode::Escape), .. } => {
                    break 'running
                },
                Event::KeyDown { keycode: Some(Keycode::Up), .. } => {
                    dir = 1;
                },
                Event::KeyDown { keycode: Some(Keycode::Down), .. } => {
                    dir = -1;
                },
                _ => {}
            }
        }

        for _ in 0..100 {
            state.step(dir);
            dir = 0;
        }
        state.draw(&mut canvas)?;
        canvas.present();
        std::thread::sleep(Duration::new(0, 1_000_000_000u32 / 60));
    };
    Ok(())
}
