# Codex task: make the Beam Data Viewer left sidebar resizable

## Target file

Modify:

```text
beam_data_viewer.py
```

This is the Tkinter/ttk beam data viewer script. The UI is built inside:

```python
def _build_ui(self) -> None:
```

The current layout uses two fixed grid columns on `root`:

```python
root.columnconfigure(0, weight=0)
root.columnconfigure(1, weight=1)
...
left = ttk.Frame(root, padding=8)
left.grid(row=0, column=0, sticky="nsw")
...
right = ttk.Frame(root, padding=(0, 8, 8, 8))
right.grid(row=0, column=1, sticky="nsew")
```

This makes the left panel too wide/fixed and not user-resizable.

## Required change

Replace the fixed two-column `root` grid layout with a horizontal `ttk.PanedWindow` so the user can drag the divider between the left sidebar and the plot area.

The left sidebar should initially use about **25% of the window width**.

The right plot area should use about **75% of the window width**.

The divider must be draggable/resizable by the user.

## Implementation requirements

### 1. Keep the status bar at the bottom

The status label should remain outside the paned window, spanning the full window width at the bottom.

Use `root` with:

```python
root.rowconfigure(0, weight=1)
root.rowconfigure(1, weight=0)
root.columnconfigure(0, weight=1)
```

The paned window should be placed at row 0, column 0.

The status label should be placed at row 1, column 0.

### 2. Add a horizontal paned window

Inside `_build_ui`, create:

```python
main_pane = ttk.PanedWindow(root, orient="horizontal")
main_pane.grid(row=0, column=0, sticky="nsew")
```

Then create the left and right frames as children of `main_pane`, not children of `root`:

```python
left = ttk.Frame(main_pane, padding=8)
right = ttk.Frame(main_pane, padding=(0, 8, 8, 8))
```

Add them to the pane:

```python
main_pane.add(left, weight=1)
main_pane.add(right, weight=3)
```

This gives an approximate 25% / 75% initial split.

### 3. Store the pane as an instance variable

Store it for future access:

```python
self.main_pane = main_pane
```

### 4. Set the initial sash position after Tk lays out the window

Because the window width may not be final during `_build_ui`, schedule the initial sash position using `after_idle` or `after`.

Add a helper method to the `BeamDataViewer` class:

```python
def _set_initial_pane_position(self) -> None:
    """Set left sidebar to ~25% of current window width on first layout."""
    try:
        width = self.root.winfo_width()
        if width <= 1:
            width = 1250
        self.main_pane.sashpos(0, int(width * 0.25))
    except Exception:
        pass
```

Then call it at the end of `_build_ui`:

```python
self.root.after_idle(self._set_initial_pane_position)
```

### 5. Remove old two-column root grid assumptions

Remove or replace these lines:

```python
root.columnconfigure(0, weight=0)
root.columnconfigure(1, weight=1)
root.rowconfigure(0, weight=1)
```

with:

```python
root.columnconfigure(0, weight=1)
root.rowconfigure(0, weight=1)
root.rowconfigure(1, weight=0)
```

Also change the final status label grid call from:

```python
ttk.Label(root, textvariable=self.status, anchor="w").grid(row=1, column=0, columnspan=2, sticky="ew", padx=8, pady=(0, 4))
```

to:

```python
ttk.Label(root, textvariable=self.status, anchor="w").grid(row=1, column=0, sticky="ew", padx=8, pady=(0, 4))
```

### 6. Do not change parser or plotting logic

Do not change these parts unless necessary:

- MCC parser
- ASC parser
- CSV parser
- normalization logic
- default autoload paths
- visibility/toggle logic
- plotting/redraw logic

Only modify layout code.

## Expected result

When the script starts:

- The left loaded-curve/control panel occupies roughly 25% of the window width.
- The right plot notebook occupies roughly 75% of the window width.
- The user can drag the vertical divider to resize the sidebar.
- The metadata table and plot tabs still resize correctly.
- The status bar remains at the bottom across the full window.

## Suggested final structure inside `_build_ui`

The start of `_build_ui` should look like this conceptually:

```python
def _build_ui(self) -> None:
    tk = self.tk
    ttk = self.ttk

    root = self.root
    root.columnconfigure(0, weight=1)
    root.rowconfigure(0, weight=1)
    root.rowconfigure(1, weight=0)

    main_pane = ttk.PanedWindow(root, orient="horizontal")
    main_pane.grid(row=0, column=0, sticky="nsew")
    self.main_pane = main_pane

    left = ttk.Frame(main_pane, padding=8)
    left.rowconfigure(2, weight=1)

    right = ttk.Frame(main_pane, padding=(0, 8, 8, 8))
    right.rowconfigure(0, weight=1)
    right.columnconfigure(0, weight=1)

    main_pane.add(left, weight=1)
    main_pane.add(right, weight=3)

    # existing left panel widgets go here
    # existing right plot notebook widgets go here

    self.status = tk.StringVar(value="Load .mcc, .asc, and/or .csv files to begin.")
    ttk.Label(root, textvariable=self.status, anchor="w").grid(
        row=1, column=0, sticky="ew", padx=8, pady=(0, 4)
    )

    self.root.after_idle(self._set_initial_pane_position)
```

## Validation steps

Run:

```bash
python3 -m py_compile beam_data_viewer.py
python3 beam_data_viewer.py
```

Then confirm manually:

1. The app opens without traceback.
2. The left panel starts at roughly 25% of the window width.
3. Dragging the divider resizes the left and right panes.
4. The Auto-load, 10x10 compare, show/hide toggles, and plot tabs still work.
