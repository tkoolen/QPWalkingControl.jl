function InteractBase.Widget(controller::PushApplier; max_force::Number, max_Δt::Number)
    pushbutton = button("Apply push")
    angleslider = slider(range(Float64(-π), stop=Float64(π), length=51), label="θ")
    magnitudeslider = slider(range(0.0, stop=max_force, length=51), label="magnitude")
    timeslider = slider(range(0, stop=max_Δt, length=51), label="Δt")

    let controller = controller
        Observables.on(pushbutton) do _
            controller.new_push = true
        end
        Observables.onany(angleslider, magnitudeslider) do θ, magnitude
            frame = controller.jacobian.frame
            s, c = sincos(θ)
            controller.force = FreeVector3D(frame, magnitude * c, magnitude * s, 0.0)
        end
        Observables.on(timeslider) do Δt
            controller.Δt = Δt
        end
    end

    angleslider[] = 0
    magnitudeslider[] = max_force / 2
    timeslider[] = max_Δt

    vbox(pushbutton, angleslider, magnitudeslider, timeslider)
end

struct PushRecoveryGUI
    rgui::RigidBodySim.GUI
    pushwidget
end

PushRecoveryGUI(mvis::MechanismVisualizer, pushwidget) = PushRecoveryGUI(RigidBodySim.GUI(mvis), pushwidget)

function Base.open(gui::PushRecoveryGUI, window::Blink.Window)
    Blink.title(window, "PushRecovery sim")
    rgui = gui.rgui
    Blink.body!(window, vbox(
        WebIO.render(gui.pushwidget),
        RigidBodySim.Visualization.render_default(rgui.controls),
        rgui.visualizer.visualizer.core
    ))
    wait(rgui)
    nothing
end

Base.open(gui::PushRecoveryGUI) = open(gui, Blink.Window())

RigidBodySim.CallbackSet(gui::PushRecoveryGUI) = RigidBodySim.CallbackSet(gui.rgui)
