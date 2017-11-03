use std::f64;
use std::collections::HashMap;
use std::ops::{Index, IndexMut};

use glutin;
use glutin::{ElementState, EventsLoop, VirtualKeyCode};
use glutin::WindowEvent::*;

use world::*;

pub const INPUT_UP: usize = 0;
pub const INPUT_DOWN: usize = 1;
pub const INPUT_LEFT: usize = 2;
pub const INPUT_RIGHT: usize = 3;

pub struct Input {
    pub move_forward: bool,
    pub move_backward: bool,
    pub strafe_left: bool,
    pub strafe_right: bool,
    pub action: bool,
    pub delta_x: f64,
    pub delta_y: f64,
    max_delta_len: f64,
    bindings: HashMap<VirtualKeyCode, usize>,
}


impl Input {
    pub fn new() -> Self {
        Input {
            move_forward: false,
            move_backward: false,
            strafe_left: false,
            strafe_right: false,
            action: false,
            delta_x: 0.0,
            delta_y: 0.0,
            max_delta_len: 0.0,
            bindings: HashMap::new(),
        }
    }

    fn set_mouse(&mut self, x: f64, y: f64) {
        let new_delta_x = x - SCREEN_WIDTH as f64 * 0.5;
        let new_delta_y = y - SCREEN_HEIGHT as f64 * 0.5;
        let dot = new_delta_x * new_delta_x + new_delta_y * new_delta_y;
        if dot > self.max_delta_len {
            self.delta_x = new_delta_x;
            self.delta_y = new_delta_y;
            self.max_delta_len = dot
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
        self.max_delta_len = 0.0;
//        let bindings = &self.bindings;
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
            _ => &self.action,
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
            _ => &mut self.action,
        }
    }
}

/*
/// MenuInput describes buttons that could be pressed to perform menu actions,
/// plus the location of the mouse.
pub struct MenuInput {
    pub option_up: bool,
    pub option_down: bool,
    pub option_left: bool,
    pub option_right: bool,
    pub clicked: bool,
    pub mouse_x: f64,
    pub mouse_y: f64,
}

impl MenuInput {
    fn new() -> Self {
        MenuInput {
            option_up: false,
            option_down: false,
            option_left: false,
            option_right: false,
            clicked: false,
            mouse_x: 0.0,
            mouse_y: 0.0,
        }
    }
}

impl Index<usize> for MenuInput {
    type Output = bool;

    fn index<'a>(&'a self, index: usize) -> &'a bool {
        match index {
            INPUT_UP => &self.option_up,
            INPUT_DOWN => &self.option_down,
            INPUT_LEFT => &self.option_left,
            INPUT_RIGHT => &self.option_right,
            _ => &self.clicked,
        }
    }
}

impl IndexMut<usize> for MenuInput {
    fn index_mut<'a>(&'a mut self, index: usize) -> &'a mut bool {
        match index {
            INPUT_UP => &mut self.option_up,
            INPUT_DOWN => &mut self.option_down,
            INPUT_LEFT => &mut self.option_left,
            INPUT_RIGHT => &mut self.option_right,
            _ => &mut self.clicked,
        }
    }
}

impl InputScenario for MenuInput {
    fn enter_frame(&mut self) {}

    fn set_mouse(&mut self, x: f64, y: f64) {
        self.mouse_x = x;
        self.mouse_y = y;
    }
}

pub struct InputState {
    // Maybe make an enum, use that to store both
    pub game_state: GameInput,
    pub menu_state: MenuInput,
    bindings: HashMap<VirtualKeyCode, usize>,
    curr_state: GameState,
}

impl InputState {
    pub fn new() -> Self {
        InputState {
            game_state: GameInput::new(),
            menu_state: MenuInput::new(),
            bindings: HashMap::new(),
            curr_state: GameState::Running,
        }
    }

    pub fn bind_key(&mut self, key: VirtualKeyCode, action: usize) {
        self.bindings.insert(key, action);
    }

    pub fn gather(&mut self, events_loop: &mut EventsLoop) -> GameState {
        self.game_state.enter_frame();
        self.menu_state.enter_frame();
        let curr_input: &mut InputScenario<Output=bool> = match self.curr_state {
            GameState::MainMenu | GameState::Paused | GameState::Exiting =>
                &mut self.menu_state,
            _ =>
                &mut self.game_state,
        };
        let curr_state = &mut self.curr_state;
        let bindings = &self.bindings;
        events_loop.poll_events(|event|{
            if let glutin::Event::WindowEvent{event, ..} = event {
                match event {
                    Closed |
                    KeyboardInput{ input: glutin::KeyboardInput{ virtual_keycode: Some(VirtualKeyCode::Escape), .. }, .. } => 
                    //KeyboardInput(_, _, Some(VirtualKeyCode::Escape), _) =>
                        *curr_state = GameState::Exiting,
                    KeyboardInput{ input: glutin::KeyboardInput{ virtual_keycode: Some(key), state, .. }, .. } => {
                        if let Some(action) = bindings.get(&key) {
                            curr_input[*action] = ElementState::Pressed == state;
                        }
                    },
                    MouseMoved{ position: (x, y), .. } => curr_input.set_mouse(x, y),
                    _ => {},
                }
            }
        });
        *curr_state
    }
}
*/
