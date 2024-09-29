use std::time::Duration;
//use std::ops::{Index,IndexMut,Mul,Add,AddAssign};
use rand::Rng;
use std::iter;

use nalgebra::{SMatrix,SVector};

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

type Vector = SVector<f64, NDF>;
type Matrix = SMatrix<f64, NDF, NDF>;
type StateVector = SMatrix<f64, 2, NDF>;

struct State {
    s: StateVector,
}

impl State {
    fn sv_from_halves(x: &Vector, v: &Vector) -> StateVector {
        StateVector::from_iterator(
            x.iter().cloned().chain(v.iter().cloned())
        )
    }
    fn from_halves(x: &Vector, v: &Vector) -> State {
        State {s: Self::sv_from_halves(x, v)}
    }

    fn from_rand() -> State {
        let mut rng = rand::thread_rng();
        State::from_halves(
            &Vector::from_fn(|_, _| 8.*rng.gen::<f64>() - 0.5),
            &Vector::from_fn(|_, _| 4.*rng.gen::<f64>() - 0.5),
        )
    }

    fn deriv(p: StateVector, dir: i32) -> StateVector {
        // s_{jk} = min(j+1, k+1)*sin(phi_j-phi_k)
        // c_{jk} = min(j+1, k+1)*cos(phi_j-phi_k)
        let mut s = Matrix::zeros();
        let mut c = Matrix::zeros();
        for i in 0..NDF {
            for j in 0..i {
                s[(i,j)] = (j+1) as f64*(p[i]-p[j]).sin();
                s[(j,i)] = -s[(i,j)];
                c[(i,j)] = (j+1) as f64*(p[i]-p[j]).cos();
                c[(j,i)] = c[(i,j)];
            }
        }
        // RHS vector
        let mut y = [0.; NDF];
        for i in 0..NDF {
            y[i] = 0.5 * i as f64 * p[i].sin();
            for j in 0..NDF {
                y[i] += s[(i,j)] * sqr(p[NDF+j]);
            }
        }
        // TODO: Invert c*x=y using CG

        State::sv_from_halves(
            &Vector::from_iterator(p.row(1).iter().cloned()),
            &Vector::from_iterator([0.; NDF].iter().cloned()),
        )
    }

    fn energy(&self) -> f64 {
        0.
    }

    fn step(&mut self, dir: i32) {
        let h = 0.001;
        let k1 = Self::deriv(self.s, dir);
        let k2 = Self::deriv(self.s + 0.5*h*k1, dir);
        let k3 = Self::deriv(self.s + 0.5*h*k2, dir);
        let k4 = Self::deriv(self.s + h*k3, dir);
        self.s += h/6.0 * (k1+2.*k2+2.*k3+k4);
    }

    fn draw(&self, canvas: &mut Canvas<Window>) -> Result<(), String> {
        println!("{}", self.energy());
        let mut pos = (512., 384.);
        let mut nextpos: (f64, f64);
        for i in (0..NDF).rev() {
            let (s, c) = self.s[i].sin_cos();
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
