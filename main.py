# main.py
from gui import DeepSkyApp
import ttkbootstrap as tb
from horizon import load_config

if __name__ == "__main__":
    cfg   = load_config()
    theme = cfg.get("theme", "darkly")

    # 1) Create the style
    style = tb.Style(theme=theme)
    # 2) Globally set every ttk widget to Helvetica 12
    style.configure('.', font=('Helvetica', 12))

    # 3) Grab the underlying Tk root and also set the menu font
    root = style.master
    root.option_add('*Menu.font', 'Helvetica 12')

    # 4) Launch
    app = DeepSkyApp(root)
    root.mainloop()
