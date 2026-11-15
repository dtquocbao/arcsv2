# Run arcsv2 via Compose (project name arcsv2). Repo is /work, data is /data.
# Quote-safe for a checkout path that contains spaces or an apostrophe.
$ErrorActionPreference = "Stop"
$Root = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot "..")).Path
$ComposeFile = Join-Path $Root "compose.yaml"

function Invoke-Compose {
    param([string[]]$ComposeArgs)
    & docker @("compose", "-f", $ComposeFile, "--project-directory", $Root) @ComposeArgs
    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
}

function Ensure-Data {
    # Named volumes are created root-owned; mambauser must be able to write /data.
    Invoke-Compose @("--profile", "init", "run", "--rm", "init-data") | Out-Null
}

function ConvertTo-ContainerPath {
    param([string]$PathArg)
    if ([string]::IsNullOrWhiteSpace($PathArg)) { return $PathArg }
    if ($PathArg.StartsWith("/")) { return $PathArg }
    try {
        $full = [System.IO.Path]::GetFullPath((Join-Path (Get-Location).Path $PathArg))
    } catch {
        return $PathArg
    }
    $rootFull = [System.IO.Path]::GetFullPath($Root)
    if ($full.StartsWith($rootFull, [StringComparison]::OrdinalIgnoreCase)) {
        $rel = $full.Substring($rootFull.Length).TrimStart("\", "/")
        if ([string]::IsNullOrEmpty($rel)) { return "/work" }
        return "/work/" + ($rel -replace "\\", "/")
    }
    return $PathArg
}

function Get-DashboardArgs {
    param([System.Collections.Generic.List[string]]$AfterDashboard)
    $out = [System.Collections.Generic.List[string]]::new()
    [void]$out.Add("dashboard")
    $hasHost = $false
    $hasResults = $false
    for ($i = 0; $i -lt $AfterDashboard.Count; $i++) {
        $a = $AfterDashboard[$i]
        if ($a -eq "--host") { $hasHost = $true }
        if ($a -eq "--results") {
            $hasResults = $true
            [void]$out.Add($a)
            if ($i + 1 -lt $AfterDashboard.Count) {
                $i++
                [void]$out.Add((ConvertTo-ContainerPath $AfterDashboard[$i]))
            }
            continue
        }
        [void]$out.Add($a)
    }
    if (-not $hasResults) {
        [void]$out.Add("--results")
        # Fetch/Snakemake write under the named volume, not the repo bind.
        [void]$out.Add("/data")
    }
    if (-not $hasHost) {
        [void]$out.Add("--host")
        [void]$out.Add("0.0.0.0")
    }
    return $out.ToArray()
}

$Rest = [System.Collections.Generic.List[string]]::new()
foreach ($a in $args) { [void]$Rest.Add([string]$a) }

if ($Rest.Count -ge 1 -and $Rest[0] -eq "build") {
    Invoke-Compose @("build", "cli")
    exit 0
}
if ($Rest.Count -ge 1 -and $Rest[0] -eq "dashboard") {
    Ensure-Data
    $after = [System.Collections.Generic.List[string]]::new()
    if ($Rest.Count -gt 1) {
        foreach ($x in $Rest.GetRange(1, $Rest.Count - 1)) { [void]$after.Add($x) }
    }
    $dashArgs = Get-DashboardArgs $after
    $cmd = [System.Collections.Generic.List[string]]::new()
    foreach ($x in @(
            "--profile", "dashboard", "run", "--rm", "--service-ports", "dashboard"
        )) { [void]$cmd.Add($x) }
    foreach ($x in $dashArgs) { [void]$cmd.Add($x) }
    Invoke-Compose $cmd.ToArray()
    exit 0
}
if ($Rest.Count -ge 1 -and $Rest[0] -eq "shell") {
    Ensure-Data
    Invoke-Compose @("--profile", "cli", "run", "--rm", "--entrypoint", "bash", "cli")
    exit 0
}
if ($Rest.Count -ge 1 -and $Rest[0] -eq "fetch") {
    Ensure-Data
    $cmd = [System.Collections.Generic.List[string]]::new()
    foreach ($x in @(
            "--profile", "cli", "run", "--rm", "--entrypoint", "python", "cli",
            "scripts/fetch_data.py", "--auto"
        )) { [void]$cmd.Add($x) }
    if ($Rest.Count -gt 1) {
        foreach ($x in $Rest.GetRange(1, $Rest.Count - 1)) { [void]$cmd.Add($x) }
    }
    Invoke-Compose $cmd.ToArray()
    exit 0
}
if ($Rest.Count -ge 2 -and $Rest[0] -eq "--entrypoint") {
    Ensure-Data
    $ep = $Rest[1]
    $cmd = [System.Collections.Generic.List[string]]::new()
    foreach ($x in @("--profile", "cli", "run", "--rm", "--entrypoint", $ep, "cli")) { [void]$cmd.Add($x) }
    if ($Rest.Count -gt 2) {
        foreach ($x in $Rest.GetRange(2, $Rest.Count - 2)) { [void]$cmd.Add($x) }
    }
    Invoke-Compose $cmd.ToArray()
    exit 0
}
Ensure-Data
$cmd = [System.Collections.Generic.List[string]]::new()
foreach ($x in @("--profile", "cli", "run", "--rm", "cli")) { [void]$cmd.Add($x) }
foreach ($x in $Rest) { [void]$cmd.Add($x) }
Invoke-Compose $cmd.ToArray()
exit 0
