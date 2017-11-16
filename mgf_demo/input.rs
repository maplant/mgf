use std::f64;
use std::collections::HashMap;
use std::ops::{Index, IndexMut};

use glutin;
use glutin::{ElementState, EventsLoop, VirtualKeyCode};
use glutin::WindowEvent::*;

pub const INPUT_UP: usize = 0;
pub const INPUT_DOWN: usize = 1;
pub const INPUT_LEFT: usize = 2;
pub const INPUT_RIGHT: usize = 3;

enum MouseState {
    Released,
    JustPressed,
    Held{ prev_x: f64, prev_y: f64 },
}

pub struct Input {
    pub move_forward: bool,
    pub move_backward: bool,
    pub strafe_left: bool,
    pub strafe_right: bool,
    mouse: MouseState,
    pub delta_x: f64,
    pub delta_y: f64,
    bindings: HashMap<VirtualKeyCode, usize>,
}


impl Input {
    pub fn new() -> Self {
        Input {
            move_forward: false,
            move_backward: false,
            strafe_left: false,
            strafe_right: false,
            mouse: MouseState::Released,
            delta_x: 0.0,
            delta_y: 0.0,
            bindings: HashMap::new(),
        }
    }

    fn set_mouse(&mut self, x: f64, y: f64) {
        match self.mouse {
            MouseState::Released => {
                self.delta_x = 0.0;
                self.delta_y = 0.0;
            },
            MouseState::JustPressed => {
                self.mouse = MouseState::Held {
                    prev_x: x,
                    prev_y: y
                };
            },
            MouseState::Held{ prev_x, prev_y } => {
                self.delta_x += x - prev_x;
                self.delta_y += y - prev_y;
                self.mouse = MouseState::Held {
                    prev_x: x,
                    prev_y: y
                };
            },
        }
    }

    pub fn bind_key(&mut self, key: VirtualKeyCode, action: usize) {
        self.bindings.insert(key, action);
    }

    fn get_binding(&self, key: &VirtualKeyCode) -> Option<usize> {
        if let Some(action) = self.bindings.get(&key) {
            Some(*action)
        } else {
            None
        }
    }

    pub fn gather(&mut self, events_loop: &mut EventsLoop) -> bool {
        self.delta_x = 0.0;
        self.delta_y = 0.0;
        let mut continue_game = true;
        events_loop.poll_events(|event|{
            if let glutin::Event::WindowEvent{event, ..} = event {
                match event {
                    Closed |
                    KeyboardInput{ input: glutin::KeyboardInput{ virtual_keycode: Some(VirtualKeyCode::Escape), .. }, .. } =>
                        continue_game = false,
                    KeyboardInput{ input: glutin::KeyboardInput{ virtual_keycode: Some(key), state, .. }, .. } => {
                        if let Some(action) = self.get_binding(&key) {
                           self[action] = ElementState::Pressed == state;
                        }
                    },
                    MouseInput{ state: glutin::ElementState::Pressed, button: glutin::MouseButton::Left, .. } => {
                        if let MouseState::Released = self.mouse {
                            self.mouse = MouseState::JustPressed;
                        }
                    },
                    MouseInput{ state: glutin::ElementState::Released, button: glutin::MouseButton::Left, .. } => {
                        self.mouse = MouseState::Released;
                    },
                    MouseMoved{ position: (x, y), .. } => self.set_mouse(x, y),
                    _ => {},
                }
            }
        });
        continue_game
    }
}

impl Index<usize> for Input {
    type Output = bool;

    fn index<'a>(&'a self, index: usize) -> &'a bool {
        match index {
            INPUT_UP => &self.move_forward,
            INPUT_DOWN => &self.move_backward,
            INPUT_LEFT => &self.strafe_left,
            INPUT_RIGHT => &self.strafe_right,
            _ => unimplemented!(),
        }
    }
}

impl IndexMut<usize> for Input {
    fn index_mut<'a>(&'a mut self, index: usize) -> &'a mut bool {
        match index {
            INPUT_UP => &mut self.move_forward,
            INPUT_DOWN => &mut self.move_backward,
            INPUT_LEFT => &mut self.strafe_left,
            INPUT_RIGHT => &mut self.strafe_right,
            _ => unimplemented!(),
        }
    }
}
