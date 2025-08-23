# Connecting with NX

## Connecting with NX

Install [x2goclient](https://wiki.x2go.org/doku.php/doc:installation:x2goclient). It may be available as a package for your Linux/Mac distribution, or you can download it (also for Windows) from the provided link.

To set up a new session:

*   Click on `Session`, then `New session…`.
*   Optionally, change the `Session name`.
*   Enter `10.195.17.215` in the `Host` field.
*   Fill in the `Login` field with your username.
*   Set the location of your private key in the `Use RSA/DSA key for ssh connection` field.
*   In the `Session type` box, choose either:
    *   `XFCE` for a graphical desktop environment.
    *   `Single application`, then `Terminal` from the right drop-down menu for a terminal session.
*   Click `Ok`. The session should appear in the right panel. Click on it to connect.

![X2GO](images/Nx.png)

**Important**: If you are outside the university’s network, you must first connect to the VPN. Refer to the [Linux/Mac SSH connection tutorial](hpc_connect_ssh_linux_mac.md) or [Windows SSH connection tutorial](hpc_connect_ssh_windows.md) for VPN instructions.

