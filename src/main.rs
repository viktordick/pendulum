use std::time::Duration;
use std::ops::{Index,IndexMut,Mul,Add,AddAssign};
use rand::Rng;

use sdl2::event::Event;
use sdl2::keyboard::Keycode;
use sdl2::pixels::Color;
use sdl2::render::Canvas;
use sdl2::video::Window;
use sdl2::gfx::primitives::DrawRenderer;

fn sqr(x: f64) -> f64 {
    x*x
}

#[derive(Clone, Copy)]
struct Vector {
    f: [f64; 4],
}

impl Vector {
    fn new(x1: f64, x2: f64, x3: f64, x4: f64) -> Vector {
        Vector{ f: [x1, x2, x3, x4] }
    }
    fn zero() -> Vector {
        Vector{ f: [0.; 4] }
    }
}
impl Index<usize> for Vector {
    type Output = f64;
    fn index(&self, idx: usize) -> &f64 {
        &self.f[idx]
    }
}
impl IndexMut<usize> for Vector {
    fn index_mut(&mut self, idx: usize) -> &mut f64 {
        &mut self.f[idx]
    }
}
impl Mul<f64> for Vector {
    type Output = Self;
    fn mul(self, rhs: f64) -> Self {
        let mut result = Self::zero();
        for i in 0..4 {
            result[i] = self[i] * rhs;
        }
        result
    }
}
impl Mul<Vector> for f64 {
    type Output = Vector;
    fn mul(self, rhs: Vector) -> Vector {
        rhs*self
    }
}
impl Add for Vector {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        let mut result = Self::zero();
        for i in 0..4 { result[i] = self[i] + rhs[i] };
        result
    }
}
impl AddAssign for Vector {
    fn add_assign(&mut self, rhs: Self) {
        for i in 0..4 { self[i] += rhs[i] };
    }
}

struct State {
    p: Vector,
}

impl State {
    fn new() -> State {
        let mut rng = rand::thread_rng();
        State {
            p: Vector::new(
               8.*(rng.gen::<f64>()-0.5),
               8.*(rng.gen::<f64>()-0.5),
               4.*rng.gen::<f64>()-0.5,
               4.*rng.gen::<f64>()-0.5,
           ),
        }
    }
    fn deriv(p: Vector, dir: i32) -> Vector {
        let (s, c) = (p[1]-p[0]).sin_cos();
        let s1 = p[0].sin();
        let s2 = p[1].sin();
        let x = 2.*s1-s*sqr(p[3]);
        let y = s*sqr(p[2])+s2;
        let inv = -1./(1.+sqr(s));
        Vector::new(
            p[2],
            p[3],
            inv*(x-c*y) + dir as f64 * p[2],
            inv*(2.*y-c*x) + dir as f64 * p[3],
        )

    }
    fn energy(&self) -> f64 {
        sqr(self.p[2]) + sqr(self.p[3])/2.
            + (self.p[1]-self.p[0]).cos()*self.p[2]*self.p[3]
            - 2.*self.p[0].cos()
            - self.p[1].cos()
            + 3.
    }
    fn step(&mut self, dir: i32) {
        let h = 0.01;
        let k1 = Self::deriv(self.p, dir);
        let k2 = Self::deriv(self.p + 0.5*h*k1, dir);
        let k3 = Self::deriv(self.p + 0.5*h*k2, dir);
        let k4 = Self::deriv(self.p + h*k3, dir);
        self.p += h/6.0 * (k1+2.*k2+2.*k3+k4);
    }
    fn draw(&self, canvas: &mut Canvas<Window>) -> Result<(), String> {
        println!("{}", self.energy());
        let mut c = [(512.,384.), self.p[0].sin_cos(), self.p[1].sin_cos()];
        for i in 0..2 {
            c[i+1] = (c[i].0 + 200. * c[i+1].0, c[i].1 + 200. * c[i+1].1);
            canvas.thick_line(
                c[i].0 as i16, c[i].1 as i16,
                c[i+1].0 as i16, c[i+1].1 as i16,
                10, Color::RGB(0,0,0),
            )?;
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
        .present_vsync()
        .accelerated()
        .build()
        .map_err(|e| e.to_string())?;
    let mut event_pump = sdl_context.event_pump()?;

    let mut state = State::new();
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

        for _ in 0..10 {
            state.step(dir);
            dir = 0;
        }
        state.draw(&mut canvas)?;
        canvas.present();
        std::thread::sleep(Duration::new(0, 1_000_000_000u32 / 60));
    };
    Ok(())
}
